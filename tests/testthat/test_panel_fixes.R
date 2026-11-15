library(levi)

# Regression tests for the defects found by the expert panel (2026-11).

hub_net <- system.file("extdata", "hub_network.dat", package = "levi")
genes <- c("HUB", paste0("N", 1:8))

# --- leviSingleCellInteractionTFCE ------------------------------------------

sc_data <- function(seed, main = 1, interaction = 1, extra_genes = 40,
                    n_donors = 5) {
    set.seed(seed)
    all_genes <- c(genes, if (extra_genes) paste0("x", seq_len(extra_genes)))
    cells <- expand.grid(cell = 1:4, type = c("A", "B"),
                         cond = c("ctrl", "stim"), donor = paste0("d", seq_len(n_donors)))
    mu <- outer(rep(30, length(all_genes)), rep(1, nrow(cells)))
    stim <- cells$cond == "stim"
    mu[seq_along(genes), stim] <- mu[seq_along(genes), stim] * main
    hit <- stim & cells$type == "B"
    mu[seq_along(genes), hit] <- mu[seq_along(genes), hit] * interaction
    counts <- matrix(rpois(length(mu), mu), nrow(mu),
                     dimnames = list(all_genes, NULL))
    list(counts = counts, donor = as.character(cells$donor),
         type = as.character(cells$type), cond = as.character(cells$cond))
}

run_interaction <- function(d, ...) {
    suppressWarnings(leviSingleCellInteractionTFCE(d$counts, d$donor, d$type,
        d$cond, hub_net, reference_cell_type = "A", min_cells = 1,
        normalize = "none", permutation_method = "exact", ...))
}

test_that("interaction TFCE is unaffected by a main effect shared by all types", {
    skip_if_not_installed("limma")
    base <- run_interaction(sc_data(11))
    # The same counts with the network genes of every stim cell multiplied by
    # 4 in both cell types: the per-donor interaction contrasts are unchanged
    # up to the pseudocount and library size, so the T statistics barely move
    # and nothing becomes significant.
    d <- sc_data(11)
    stim <- d$cond == "stim"
    d$counts[genes, stim] <- d$counts[genes, stim] * 4L
    shifted <- run_interaction(d)
    expect_equal(shifted$results$B$T, base$results$B$T, tolerance = 0.2)
    expect_false(any(shifted$results$B$GlobalSignificant))
})

test_that("interaction TFCE aligns cell types that keep different genes", {
    skip_if_not_installed("limma")
    d <- sc_data(12)
    # Extra genes expressed in only one cell type are filtered from the other
    # type, so the per-type matrices have different rows.
    d$counts[paste0("x", 1:20), d$type == "A"] <- 0L
    d$counts[paste0("x", 21:40), d$type == "B"] <- 0L
    out <- run_interaction(d)
    expect_equal(out$results$B$Gene, genes)
    expect_true(all(is.finite(out$results$B$T)))
    expect_equal(length(out$null_global), 2^5 - 1)
})

test_that("interaction TFCE detects a type-specific condition effect", {
    skip_if_not_installed("limma")
    out <- run_interaction(sc_data(13, interaction = 4))
    expect_true(all(out$results$B$T > 0))
    expect_equal(min(out$results$B$PGlobal), 2 / 2^5)
})

# --- leviGraphMoran: conditional local permutation ---------------------------

lattice_net <- function() {
    g <- igraph::make_lattice(c(6, 6))
    el <- igraph::as_edgelist(g)
    list(nodes = data.frame(ID = paste0("g", 1:36), X = rep(1:6, 6),
                            Y = rep(1:6, each = 6)),
         edges = data.frame(V1 = paste0("g", el[, 1]),
                            V2 = paste0("g", el[, 2])))
}

test_that("local Moran keeps the node's own score fixed", {
    net <- lattice_net()
    set.seed(5)
    x <- stats::setNames(rnorm(36), net$nodes$ID)
    out <- suppressWarnings(leviGraphMoran(x, net$nodes, net$edges, fileTypeInput = "stg",
                          n_perm = 199, seed = 1), classes = "levi_resolution_warning")
    # Under conditional permutation the local I of a node is its own z times
    # a sum of degree-many other z values, so the null mean is
    # -k_i z_i^2 / (n - 1); with 199 draws the observed p-values of iid
    # scores should be roughly uniform.
    expect_equal(nrow(out$local), 36)
    expect_gt(mean(out$local$P), 0.3)
    # An isolated outlier with iid neighbours is not a cluster.
    x["g15"] <- 5
    out <- suppressWarnings(leviGraphMoran(x, net$nodes, net$edges, fileTypeInput = "stg",
                          n_perm = 199, seed = 1), classes = "levi_resolution_warning")
    expect_gt(out$local$P[out$local$Gene == "g15"], 0.05)
})

test_that("local Moran flags a planted cluster", {
    net <- lattice_net()
    set.seed(6)
    x <- stats::setNames(rnorm(36, sd = 0.3), net$nodes$ID)
    x[c("g14", "g15", "g16", "g9", "g21")] <- 3
    out <- suppressWarnings(leviGraphMoran(x, net$nodes, net$edges, fileTypeInput = "stg",
                          n_perm = 499, seed = 2), classes = "levi_resolution_warning")
    expect_lt(out$local$P[out$local$Gene == "g15"], 0.05)
    expect_lt(out$global$P, 0.05)
})

# --- Input validation --------------------------------------------------------

test_that("graph tests reject unnamed or unmatched scores", {
    net <- lattice_net()
    expect_error(suppressWarnings(leviGraphMoran(rnorm(36), net$nodes, net$edges,
                                fileTypeInput = "stg", n_perm = 9), classes = "levi_resolution_warning"),
                 "named by gene")
    expect_error(suppressWarnings(leviGraphMoran(c(zz = 1, yy = 2), net$nodes, net$edges,
                                fileTypeInput = "stg", n_perm = 9), classes = "levi_resolution_warning"),
                 "finite score")
})

test_that("n_perm and fileTypeInput are validated", {
    net <- lattice_net()
    x <- stats::setNames(rnorm(36), net$nodes$ID)
    for (bad in list(0, -5, 2.5, "9"))
        expect_error(suppressWarnings(leviGraphMoran(x, net$nodes, net$edges,
                                    fileTypeInput = "stg", n_perm = bad), classes = "levi_resolution_warning"),
                     "n_perm")
    expect_error(suppressWarnings(leviGraphMoran(x, net$nodes, net$edges,
                                fileTypeInput = "xyz", n_perm = 9), classes = "levi_resolution_warning"),
                 "Unknown fileTypeInput")
    expect_silent(suppressWarnings(leviGraphMoran(x, net$nodes, net$edges,
                                 fileTypeInput = "STG", n_perm = 9), classes = "levi_resolution_warning"))
})

test_that("reciprocal edges are counted once", {
    nodes <- data.frame(ID = c("A", "B", "C"), X = 1:3, Y = 1:3)
    edges <- data.frame(V1 = c("A", "B", "B"), V2 = c("B", "A", "C"))
    ni <- levi:::.networkIndex(nodes, edges, "stg")
    expect_equal(nrow(ni$edges), 2L)
})

# --- .maxPvalue --------------------------------------------------------------

test_that(".maxPvalue matches the direct count, ties included", {
    set.seed(1)
    null <- round(rnorm(500), 2)
    obs <- c(round(rnorm(100), 2), null[1:20], -10, 10)
    direct <- (1 + vapply(obs, function(o)
        sum(null >= o - 1e-12 * max(1, abs(o))), numeric(1))) / 501
    expect_equal(levi:::.maxPvalue(null, obs), direct)
})

# --- Adapters ----------------------------------------------------------------

test_that("leviFromLimma refuses the intercept and ambiguous fits", {
    skip_if_not_installed("limma")
    set.seed(2)
    y <- matrix(rnorm(9 * 6), 9, dimnames = list(genes, NULL))
    y[1:3, 4:6] <- y[1:3, 4:6] + 6
    group <- factor(rep(c("C", "T"), each = 3))
    fit <- limma::eBayes(limma::lmFit(y, stats::model.matrix(~ group)))
    expect_error(leviFromLimma(fit), "2 coefficients")
    expect_error(leviFromLimma(fit, coef = 1), "intercept")
    expect_error(leviFromLimma(fit, coef = "nope"), "not in the fit")
    out <- leviFromLimma(fit, coef = "groupT")
    expect_equal(out$GeneID, genes)
    expect_true(all(out$logFC[1:3] > 3))
})

test_that("leviFromSE accepts sparse and DelayedArray assays", {
    skip_if_not_installed("SummarizedExperiment")
    dense <- matrix(c(0, 2, 4, 0, 6, 8, 1, 0, 3, 5, 0, 7), 3,
                    dimnames = list(c("g1", "g2", "g3"), paste0("s", 1:4)))
    cd <- data.frame(condition = c("T", "T", "C", "C"),
                     row.names = colnames(dense))
    expected <- suppressWarnings(leviFromSE(SummarizedExperiment::SummarizedExperiment(
        list(counts = dense), colData = cd), condition_col = "condition",
        test_level = "T", ctrl_level = "C"))
    sparse <- Matrix::Matrix(dense, sparse = TRUE)
    se <- SummarizedExperiment::SummarizedExperiment(list(counts = sparse),
                                                     colData = cd)
    expect_equal(suppressWarnings(leviFromSE(se, condition_col = "condition", test_level = "T",
                            ctrl_level = "C")), expected)
    skip_if_not_installed("DelayedArray")
    se <- SummarizedExperiment::SummarizedExperiment(
        list(counts = DelayedArray::DelayedArray(dense)), colData = cd)
    expect_equal(suppressWarnings(leviFromSE(se, condition_col = "condition", test_level = "T",
                            ctrl_level = "C")), expected)
})

# --- Tests that fail on the code they were written against -------------------
# The tests above pass on both the old and the new code for some defects,
# because with five donors the exact p-value floor (2/32 = 0.0625) already
# exceeds 0.05. The ones below measure a rejection rate over replicates and
# are bounded at least three binomial standard errors from the old value.

test_that("interaction TFCE has a valid null under a large main effect", {
    skip_if_not_installed("limma")
    # Six donors: floor 2/2^6 = 0.031 < 0.05, so rejection is attainable. The
    # pre-1.99.5 statistic rejected in about 85% of these data sets because
    # the main effect leaked into the interaction.
    reject <- vapply(1:20, function(i) {
        d <- sc_data(100 + i, main = 4, interaction = 1, extra_genes = 10,
                     n_donors = 6)
        out <- run_interaction(d)
        any(out$results$B$GlobalSignificant)
    }, logical(1))
    expect_lte(mean(reject), 0.15)
})

test_that("interaction TFCE has power with six donors", {
    skip_if_not_installed("limma")
    reject <- vapply(1:6, function(i) {
        d <- sc_data(200 + i, main = 1, interaction = 4, extra_genes = 10,
                     n_donors = 6)
        any(run_interaction(d)$results$B$GlobalSignificant)
    }, logical(1))
    expect_gte(mean(reject), 0.8)
})

test_that("local Moran does not call an isolated outlier a cluster", {
    # Total permutation (pre-1.99.5) rejected an outlying node with iid
    # neighbours in 52% of replicates; conditional permutation in 8%.
    net <- lattice_net()
    hit <- vapply(1:60, function(i) {
        set.seed(1000 + i)
        x <- stats::setNames(rnorm(36), net$nodes$ID)
        x["g15"] <- 5
        out <- suppressWarnings(leviGraphMoran(x, net$nodes, net$edges, fileTypeInput = "stg",
                              n_perm = 99, seed = i), classes = "levi_resolution_warning")
        out$local$P[out$local$Gene == "g15"] <= 0.05
    }, logical(1))
    expect_lt(mean(hit), 0.2)
})

test_that("leviBulkGraphInference refits dispersions by default", {
    expect_true(formals(leviBulkGraphInference)$refit_dispersions)
    skip_if_not_installed("DESeq2")
    set.seed(3)
    genes <- c("HUB", paste0("N", 1:8))
    counts <- matrix(rpois(length(genes) * 6, 30), length(genes), 6,
        dimnames = list(genes, paste0("s", 1:6)))
    grp <- c("C", "C", "C", "T", "T", "T")
    net <- system.file("extdata", "hub_network.dat", package = "levi")
    fixed <- NULL
    expect_warning(
        suppressMessages(fixed <- leviBulkGraphInference(counts, grp, net,
            test = "T", control = "C", threshold = 1, n_perm = 3, seed = 1,
            permutation_method = "exact", refit_dispersions = FALSE)),
        "no level guarantee")
    expect_false(fixed$refit_dispersions)
    # 20 sample-label arrangements: the exact floor is 2/20 = 0.1, the same
    # value as in the existing bulk test, and does not depend on the mode.
    expect_equal(fixed$p_floor, 0.1)
})

test_that("bulk inference gives identical results with the same seed", {
    skip_if_not_installed("DESeq2")
    skip_on_os("windows")
    set.seed(4)
    genes <- c("HUB", paste0("N", 1:8))
    counts <- matrix(rpois(length(genes) * 6, 30), length(genes), 6,
        dimnames = list(genes, paste0("s", 1:6)))
    grp <- c("C", "C", "C", "T", "T", "T")
    net <- system.file("extdata", "hub_network.dat", package = "levi")
    run <- function(bp) suppressWarnings(suppressMessages(
        leviBulkGraphInference(counts, grp, net, test = "T", control = "C",
            threshold = 1, n_perm = 5, seed = 7, BPPARAM = bp)))
    a <- run(BiocParallel::SerialParam())
    b <- run(BiocParallel::MulticoreParam(2))
    expect_identical(a$regions$summary$PSpatial, b$regions$summary$PSpatial)
})

# --- Input checks ------------------------------------------------------------

test_that("expression with no gene in common with the network is an error", {
    skip_if_not_installed("limma")
    x <- matrix(rnorm(60), 10, 6, dimnames = list(paste0("zz", 1:10), NULL))
    expect_error(
        leviGraphTFCEInference(x, rep(c("C", "T"), 3), hub_net, n_perm = 9,
                               test = "T", control = "C"),
        "No gene")
})

test_that("duplicated gene names and unknown groups are errors", {
    skip_if_not_installed("limma")
    x <- matrix(rnorm(54), 9, 6, dimnames = list(genes, NULL))
    dup <- x; rownames(dup)[2] <- "HUB"
    g <- rep(c("C", "T"), 3)
    expect_error(leviGraphTFCEInference(dup, g, hub_net, n_perm = 9,
                                        test = "T", control = "C"),
                 "Duplicated")
    expect_error(leviGraphTFCEInference(x, g, hub_net, n_perm = 9,
                                        test = "Z", control = "C"),
                 "levels of 'groups'")
    expect_error(leviGraphTFCEInference(x, g, hub_net, n_perm = 0,
                                        test = "T", control = "C"),
                 "n_perm")
})

# --- Additive options (Phase 3) -----------------------------------------------

test_that("resolution warning has its own class and fires only below the resolution", {
    net <- lattice_net()
    x <- stats::setNames(rnorm(36), net$nodes$ID)
    expect_warning(
        leviGraphMoran(x, net$nodes, net$edges, fileTypeInput = "stg",
                       n_perm = 99, seed = 1),
        class = "levi_resolution_warning")
    # 36 nodes, alpha 0.05: one node qualifies once 1 / (n_perm + 1) <= 0.05 / 36.
    expect_no_warning(
        leviGraphMoran(x, net$nodes, net$edges, fileTypeInput = "stg",
                       n_perm = 799, seed = 1))
})

test_that("Getis-Ord reports a conditional P-value that keeps the node's own score", {
    net <- lattice_net()
    # Rejection rate among extreme nodes under an iid null: the unconditional
    # P moves the node's own value as well and rejected about 11%, the
    # conditional one about 5%.
    big_total <- big_cond <- logical(0)
    for (i in 1:40) {
        set.seed(300 + i)
        x <- stats::setNames(rnorm(36), net$nodes$ID)
        out <- suppressWarnings(leviGraphGetisOrd(x, net$nodes, net$edges,
            fileTypeInput = "stg", n_perm = 99, seed = i))
        big <- abs(as.numeric(scale(x))) > 1.5
        big_total <- c(big_total, out$PTwoSided[big] <= 0.05)
        big_cond <- c(big_cond, out$PConditional[big] <= 0.05)
    }
    expect_true(all(c("PConditional", "PConditionalAdjusted", "Null") %in%
                    names(out)))
    expect_lt(mean(big_cond), 0.12)
    expect_gt(mean(big_total), mean(big_cond))
})

test_that("sig_level moves the flags and nothing else", {
    skip_if_not_installed("limma")
    set.seed(8)
    x <- matrix(rnorm(9 * 8), 9, 8, dimnames = list(genes, NULL))
    x[1:5, 5:8] <- x[1:5, 5:8] + 3
    g <- rep(c("C", "T"), each = 4)
    run <- function(level) suppressWarnings(leviGraphTFCEInference(x, g,
        hub_net, test = "T", control = "C", n_perm = 99, seed = 2,
        sig_level = level))
    a <- run(0.05); b <- run(0.5)
    expect_identical(a$statistic$PGlobal, b$statistic$PGlobal)
    expect_true(sum(b$statistic$GlobalSignificant) >=
                sum(a$statistic$GlobalSignificant))
    expect_error(run(1.5), "sig_level")
})

test_that("test and control from the label order warn and are stored", {
    skip_if_not_installed("limma")
    set.seed(9)
    x <- matrix(rnorm(9 * 6), 9, 6, dimnames = list(genes, NULL))
    g <- rep(c("C", "T"), 3)
    expect_warning(leviGraphTFCEInference(x, g, hub_net, n_perm = 9),
                   class = "levi_group_order_warning")
    out <- suppressWarnings(leviGraphTFCEInference(x, g, hub_net, n_perm = 9,
                                                   test = "T", control = "C"))
    expect_equal(out$test, "T")
    expect_equal(out$control, "C")
})

# --- Degree strata and degree-normalised energy ------------------------------

star_plus_net <- function() {
    # a hub of degree 12 joined to a ring of 20 nodes: very unequal degrees
    nodes <- data.frame(ID = paste0("g", 1:33), X = rnorm(33), Y = rnorm(33))
    ring <- cbind(paste0("g", 1:20), paste0("g", c(2:20, 1)))
    hub <- cbind("g21", paste0("g", 1:12))
    chain <- cbind(paste0("g", 21:32), paste0("g", 22:33))
    list(nodes = nodes, edges = data.frame(V1 = c(ring[, 1], hub[, 1], chain[, 1]),
                                           V2 = c(ring[, 2], hub[, 2], chain[, 2])))
}

test_that("Spectrum reports a degree-normalised energy, and the old columns are unchanged", {
    set.seed(1)
    net <- star_plus_net()
    x <- stats::setNames(rnorm(33), net$nodes$ID)
    out <- suppressWarnings(leviGraphSpectrum(x, net$nodes, net$edges,
        fileTypeInput = "stg", n_perm = 49, seed = 3))
    expect_true(all(c("EnergyP", "LowFrequencyP", "DegreeNormalisedEnergy",
                      "DegreeNormalisedEnergyP", "Null") %in% names(out)))
    expect_equal(out$Null, "node_label_total")
    expect_true(out$DegreeNormalisedEnergy > 0)
})

test_that("strata permute within groups and the degree warning has its own class", {
    set.seed(2)
    net <- star_plus_net()
    # scores that follow degree
    deg <- as.numeric(table(factor(c(net$edges$V1, net$edges$V2),
                                   levels = net$nodes$ID)))
    x <- stats::setNames(deg + rnorm(33, sd = 0.3), net$nodes$ID)
    expect_warning(leviGraphSpectrum(x, net$nodes, net$edges,
        fileTypeInput = "stg", n_perm = 9), class = "levi_degree_warning")
    expect_no_warning(leviGraphSpectrum(x, net$nodes, net$edges,
        fileTypeInput = "stg", n_perm = 9, strata = "degree"))
    # draws stay inside their strata
    strata <- rep(c("a", "b", "c"), length.out = 12)
    draws <- levi:::.nodeShuffles(1:12, 40, strata)
    expect_true(all(vapply(draws, function(d)
        all(vapply(split(seq_along(d), strata), function(i)
            setequal(d[i], i), logical(1))), logical(1))))
    expect_error(leviGraphSpectrum(x, net$nodes, net$edges, fileTypeInput = "stg",
                                   n_perm = 9, strata = c(1, 2)), "strata")
    out <- leviGraphSpectrum(x, net$nodes, net$edges, fileTypeInput = "stg",
                             n_perm = 49, seed = 1, strata = "degree")
    expect_equal(out$Null, "node_label_stratified")
})
