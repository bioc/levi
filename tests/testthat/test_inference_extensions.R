test_that("pathway adjustment adds a family-wise correction", {
    make_result <- function(p) list(regions = list(summary = data.frame(PSpatial = p)))
    out <- leviAdjustPathways(list(a = make_result(c(.01, .2)),
                                   b = make_result(.03)), method = "holm")
    expect_equal(out$a$regions$summary$PFamily, c(.03, .2))
    expect_equal(out$b$regions$summary$PFamily, .06)
    expect_false(out$b$regions$summary$FamilySignificant)
})

test_that("stratified node permutations do not cross strata", {
    values <- cbind(c(-2, -1, 1, 2), c(0, 0, 0, 0))
    edge <- matrix(c(1, 2), 1)
    observed <- matrix(.5, 2, 2)
    regions <- levi:::.extractLandscapeRegions(observed, 0, 1, .1, 1)
    testthat::local_mocked_bindings(landscape_m1 = function(coord, signal, ...) {
        expect_true(all(signal[1:2] < .5))
        expect_true(all(signal[3:4] > .5))
        observed
    }, .package = "levi")
    sig <- levi:::.networkSignals(values, edge, FALSE, "logfc", 1)
    levi:::.permutationPvalues(matrix(0, 5, 2), sig$signal, sig$test,
        sig$control, observed, 2, 0, 1, 1, .05, 4,
        node_values = values, edge_index = edge, signal_mode = "logfc",
        regions = regions, perm_strata = c("low", "low", "high", "high"))
})

test_that("region gene attribution ranks local non-neutral nodes", {
    x <- list(
        regions = list(cells = data.frame(Region = "over_01", X = 0, Y = 0)),
        metadata = list(node_coordinates = rbind(c(0, 0), c(1, 1)),
            node_signal = c(.8, .9), nodes = data.frame(V1 = c("A", "B")),
            grid = list(sigma = 1, increase = .1)))
    out <- leviRegionGenes(x)
    expect_equal(out$Gene[1], "A")
    expect_true(out$Contribution[1] > out$Contribution[2])
})

test_that("replicate inference permutes sample labels and returns regional p-values", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    genes <- c("HUB", paste0("N", 1:8))
    mat <- matrix(0, length(genes), 4, dimnames = list(genes, NULL))
    mat[, 3:4] <- 1
    set.seed(4)
    expect_warning(out <- leviReplicateInference(mat, c("C", "C", "T", "T"),
        test = "T", control = "C", n_perm = 3, seed = 2,
        permutation_method = "exact",
        networkCoordinatesInput = system.file("extdata", "hub_network.dat", package = "levi"),
        fileTypeInput = "dat", resolutionValueInput = 1,
        region_threshold = .02), "smallest attainable p-value")
    expect_equal(out$regions$inference$method, "sample_label")
    expect_true(all(out$regions$summary$PSpatial >= .25))
    expect_equal(out$metadata$n_perm, 5L)
    expect_true(out$metadata$permutation_exact)
})

test_that("graph cluster inference uses network-connected limma clusters", {
    skip_if_not_installed("limma")
    genes <- c("HUB", paste0("N", 1:8))
    x <- matrix(rnorm(length(genes) * 6, sd = .1), length(genes), 6,
                dimnames = list(genes, NULL))
    x[, 4:6] <- x[, 4:6] + 2
    expect_warning(out <- leviGraphClusterInference(x,
        c("C", "C", "C", "T", "T", "T"),
        system.file("extdata", "hub_network.dat", package = "levi"),
        test = "T", control = "C", threshold = 1, n_perm = 3,
        blocks = c("d1", "d2", "d3", "d1", "d2", "d3"), seed = 1,
        permutation_method = "exact"),
        "smallest attainable p-value")
    expect_true(nrow(out$regions$summary) > 0)
    expect_true(all(out$regions$summary$Nodes >= 1))
    expect_true(all(out$regions$summary$PSpatial >= .05))
    expect_true(out$exact)
})

test_that("network statistics provide TFCE, Moran, spectrum and rewiring nulls", {
    skip_if_not_installed("limma")
    net <- system.file("extdata", "hub_network.dat", package = "levi")
    genes <- c("HUB", paste0("N", 1:8))
    x <- matrix(rnorm(genes |> length() * 4), length(genes), 4,
                dimnames = list(genes, NULL))
    expect_warning(tfce <- leviGraphTFCEInference(x, c("C", "C", "T", "T"), net,
        test = "T", control = "C", n_perm = 2, seed = 1),
        "smallest attainable p-value")
    expect_equal(nrow(tfce$statistic), length(genes))
    expect_true(all(c("PGlobal", "GlobalSignificant") %in% names(tfce$statistic)))
    score <- setNames(rnorm(length(genes)), genes)
    expect_true(is.finite(suppressWarnings(
        leviGraphMoran(score, net, n_perm = 3))$global$MoranI))
    expect_true(is.finite(leviGraphSpectrum(score, net, n_perm = 3)$LaplacianEnergy))
    expect_true(is.list(leviGraphRewiringInference(score, net, n_perm = 2)))
})

test_that("single-cell TFCE uses joint donor permutations across cell types", {
    skip_if_not_installed("limma")
    genes <- c("HUB", paste0("N", 1:8)); donor <- rep(paste0("d", 1:4), each = 4)
    type <- rep(rep(c("A", "B"), each = 2), 4)
    counts <- matrix(rpois(length(genes) * length(donor), 10), length(genes),
        dimnames = list(genes, NULL))
    expect_warning(out <- leviSingleCellTFCEInference(counts, donor, type,
        c(d1 = "C", d2 = "C", d3 = "T", d4 = "T"),
        system.file("extdata", "hub_network.dat", package = "levi"), min_cells = 2,
        permutation_method = "exact", seed = 1), "smallest attainable p-value")
    expect_true(out$exact)
    expect_equal(length(out$null_global), 5L)
    expect_true(all(c("PGlobal", "GlobalSignificant") %in% names(out$results$A)))
})

test_that("single-cell TFCE respects paired donor-by-condition pseudobulks", {
    skip_if_not_installed("limma")
    genes <- c("HUB", paste0("N", 1:8))
    donor <- rep(paste0("d", 1:4), each = 4)
    condition <- rep(rep(c("ctrl", "stim"), each = 2), 4)
    type <- rep(c("A", "B"), times = 8)
    counts <- matrix(rpois(length(genes) * length(donor), 20), length(genes),
        dimnames = list(genes, NULL))
    expect_warning(out <- leviSingleCellTFCEInference(counts, donor, type,
        condition,
        system.file("extdata", "hub_network.dat", package = "levi"), min_cells = 1,
        permutation_method = "exact", seed = 1), "smallest attainable p-value")
    expect_true(out$exact)
    expect_equal(length(out$null_global), 15L)
    expect_match(out$method, "within-donor")
})

test_that("pseudobulk sums cells by donor and cell type", {
    x <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2,
                dimnames = list(c("A", "B"), NULL))
    out <- leviPseudobulk(x, donor = c("d1", "d1", "d2"),
        cell_type = c("T", "T", "B"), min_cells = 1)
    expect_equal(unname(out$counts[, "d1_T"]), c(4, 6))
    expect_equal(unname(out$counts[, "d2_B"]), c(5, 6))
})

test_that("bulk RNA-seq graph inference returns regional P-values", {
    skip_if_not_installed("DESeq2")
    set.seed(3)
    genes <- c("HUB", paste0("N", 1:8))
    counts <- matrix(rpois(length(genes) * 6, 30), length(genes), 6,
        dimnames = list(genes, paste0("s", 1:6)))
    counts[, 4:6] <- counts[, 4:6] + 30L
    expect_warning(out <- leviBulkGraphInference(counts, c("C", "C", "C", "T", "T", "T"),
        system.file("extdata", "hub_network.dat", package = "levi"),
        test = "T", control = "C", threshold = 1, n_perm = 3, seed = 1,
        permutation_method = "exact"),
        "tied pairs")
    expect_equal(out$p_floor, 0.1)
    expect_true(all(out$regions$summary$PSpatial >= .25))
})

test_that("pseudobulk expression drops unexpressed genes but keeps network nodes", {
    skip_if_not_installed("limma")
    genes <- c("HUB", paste0("N", 1:8)); donor <- rep(paste0("d", 1:4), each = 4)
    condition <- rep(rep(c("ctrl", "stim"), each = 2), 4)
    type <- rep("A", 16)
    set.seed(3)
    counts <- matrix(rpois(length(genes) * length(donor), 20), length(genes),
        dimnames = list(genes, NULL))
    counts["N8", ] <- 0                           # a network node never detected
    zeros <- matrix(0L, 2000, ncol(counts),        # a droplet-style zero block
        dimnames = list(paste0("Z", 1:2000), NULL))
    counts <- rbind(counts, zeros)
    net <- system.file("extdata", "hub_network.dat", package = "levi")
    # Fitting the 2000 all-zero rows made limma warn "eBayes unreliable" on
    # every permutation; they are filtered out now, so no such warning.
    warns <- character()
    out <- withCallingHandlers(
        leviSingleCellTFCEInference(counts, donor, type, condition, net,
            min_cells = 1, permutation_method = "exact", seed = 1),
        warning = function(w) {
            warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
        })
    expect_false(any(grepl("eBayes unreliable", warns)))
    expect_true(any(grepl("smallest attainable p-value", warns)))
    # Every network node, including the undetected one, keeps a statistic.
    expect_setequal(out$results$A$Gene, genes)
    pb <- leviPseudobulk(counts, donor, type, condition, min_cells = 1)
    keys <- paste(pb$donor, pb$condition, sep = "::")
    expr <- levi:::.pseudobulkExpression(pb, keys, "A", paired = TRUE,
        normalize = "none", keep_genes = genes)$A
    expect_true(all(genes %in% rownames(expr)))
    expect_false(any(grepl("^Z", rownames(expr))))
})

# -----------------------------------------------------------------------------
# Attainable p-value floor: 2/n under exact enumeration of a sign-symmetric
# design, 1/(n_perm + 1) under Monte Carlo; ties counted with tolerance.
# -----------------------------------------------------------------------------

test_that("3 vs 3 unblocked exact design warns at alpha = 0.05 with floor 0.10", {
    # choose(6, 3) = 20 arrangements come in complement pairs that tie the
    # maximum over both directions, so the floor is 2/20, above 0.05.
    expect_warning(
        perms <- levi:::.labelPermutations(rep(c("a", "b"), each = 3)),
        "20 distinct label arrangements.*tied pairs.*0\\.100")
    expect_true(perms$exact)
    expect_equal(perms$possible, 20)
    expect_equal(perms$floor, 0.1)
    # 4 vs 4 (70 arrangements, floor 2/70 = 0.029) can still reach 0.05.
    expect_no_warning(p4 <- levi:::.labelPermutations(rep(c("a", "b"), each = 4)))
    expect_equal(p4$floor, 2 / 70)
    # An unbalanced exact design has no complement in the enumeration, so
    # its floor stays 1/n_possible (15 arrangements: 0.067, still above 0.05).
    expect_warning(p24 <- levi:::.labelPermutations(rep(c("a", "b"), c(2, 4))),
        "15 distinct label arrangements are possible, so")
    expect_equal(p24$floor, 1 / 15)
    # Monte Carlo keeps 1/(n_perm + 1).
    expect_warning(mc <- levi:::.labelPermutations(rep(c("a", "b"), each = 10),
        n_perm = 9L, method = "monte_carlo"), "Monte Carlo.*0\\.100")
    expect_false(mc$exact)
    expect_equal(mc$floor, 0.1)
})

test_that("four paired blocks report a floor of 2/16 = 0.125", {
    groups <- rep(c("C", "T"), 4)
    blocks <- rep(paste0("d", 1:4), each = 2)
    expect_warning(perms <- levi:::.labelPermutations(groups, blocks),
        "16 distinct label arrangements.*tied pairs.*0\\.125")
    expect_equal(perms$possible, 16)
    expect_equal(perms$floor, 0.125)
    expect_length(perms$labels, 15L)
})

test_that(".maxPvalue counts a null maximum that ties up to rounding", {
    null_max <- c(1, 1 - 1e-14, 0.5)
    expect_equal(levi:::.maxPvalue(null_max, 1), 3 / 4)
    expect_equal(levi:::.maxPvalue(null_max, c(1, 0.5, 2)), c(3, 4, 1) / 4)
    expect_equal(levi:::.maxPvalue(c(1e6, 1e6 * (1 - 1e-13)), 1e6), 3 / 3)
})

test_that("exact two-group TFCE and cluster p-values never fall below 2/n", {
    skip_if_not_installed("limma")
    net <- system.file("extdata", "hub_network.dat", package = "levi")
    genes <- c("HUB", paste0("N", 1:8))
    groups <- rep(c("C", "T"), each = 3)   # 20 arrangements, floor 0.10
    for (seed in 1:5) {
        set.seed(seed)
        x <- matrix(rnorm(length(genes) * 6, sd = .2), length(genes), 6,
                    dimnames = list(genes, NULL))
        x[c("HUB", "N1", "N2"), 4:6] <- x[c("HUB", "N1", "N2"), 4:6] + 3
        tfce <- suppressWarnings(leviGraphTFCEInference(x, groups, net,
            test = "T", control = "C", permutation_method = "exact", seed = 1))
        expect_true(tfce$exact)
        expect_equal(tfce$possible_permutations, 20)
        expect_equal(tfce$p_floor, 0.1)
        # The floor applies to the joint (both-direction) p-value; the
        # per-direction POver/PUnder are compared with one direction's null,
        # where the complement is not a tie, and may reach 1/n.
        expect_true(all(tfce$statistic$PGlobal >= 0.1 - 1e-12))
        expect_true(any(tfce$statistic$PGlobal == 0.1))   # planted module is at the floor
        expect_true(all(tfce$statistic$POver >= 0.05 - 1e-12))
        cl <- suppressWarnings(leviGraphClusterInference(x, groups, net,
            test = "T", control = "C", threshold = 1,
            permutation_method = "exact", seed = 1))
        expect_equal(cl$p_floor, 0.1)
        expect_true(all(cl$regions$summary$PSpatial >= 0.1 - 1e-12))
    }
})

test_that("leviFromSeurat says its output is for visualisation, not donor-level inference", {
    markers <- data.frame(avg_log2FC = c(1.5, -0.8), pct.2 = c(0.1, 0.4),
                          row.names = c("GENE1", "GENE2"))
    expect_warning(out <- leviFromSeurat(markers), "leviPseudobulk")
    expect_equal(out$Test, c(1.5, -0.8))
})

test_that("single-cell landscape null uses the observed logfc_k", {
    skip_if_not_installed("limma")
    set.seed(1)
    genes <- c("HUB", paste0("N", 1:8))
    donor <- rep(paste0("d", 1:5), each = 4)
    condition <- rep(rep(c("ctrl", "stim"), each = 2), 5)
    cell_type <- rep(c("A", "B"), times = 10)
    counts <- matrix(rpois(length(genes) * length(donor), 20), length(genes),
                     dimnames = list(genes, NULL))
    counts["HUB", condition == "stim"] <- counts["HUB", condition == "stim"] * 4
    net <- system.file("extdata", "hub_network.dat", package = "levi")
    # The exact null contains the all-donor swap, whose landscape mirrors the
    # observed one; with one cell type its joint maximum must equal the
    # observed maximum mass, which holds only if the null uses the same k.
    for (k in c(0.5, 3)) {
        out <- suppressWarnings(leviSingleCellRegionalInference(counts, donor,
            cell_type, condition, net, min_cells = 1, threshold = 1,
            cell_types = "A", permutation_method = "exact", seed = 1,
            resolutionValueInput = 1, region_threshold = .02, logfc_k = k))
        reg <- out$landscapes$A$regions
        observed <- max(reg$summary$Mass)
        expect_true(any(abs(reg$null_max_mass$both - observed) <=
                        1e-8 * max(1, observed)))
    }
})


test_that("pseudobulk preserves labels and distinct tuples containing delimiters", {
    x <- matrix(1:8, 2, dimnames = list(c("A", "B"), NULL))
    donor <- c("d::1", "d", "d::1", "d")
    type <- c("T", "1::T", "T", "1::T")
    out <- leviPseudobulk(x, donor, type, min_cells = 1)
    expect_equal(ncol(out$counts), 2L)
    for (i in seq_along(out$donor)) {
        keep <- donor == out$donor[i] & type == out$cell_type[i]
        expect_equal(unname(out$counts[, i]), unname(rowSums(x[, keep, drop = FALSE])))
        expect_equal(unname(out$n_cells[i]), sum(keep))
    }
    expect_error(leviPseudobulk(x, donor, type, rep(NA_character_, 4)),
                 "must not contain NA")
    condition <- c("C::1", "C", "C::1", "C")
    paired <- leviPseudobulk(x, donor, type, condition, min_cells = 1)
    expect_equal(paired$condition, condition[match(paired$donor, donor)])
    expect_error(leviPseudobulk(x, c(NA, donor[-1]), type), "must not contain NA")
    collision <- leviPseudobulk(x[, 1:2], c("a_b", "a"), c("c", "b_c"), min_cells = 1)
    expect_equal(anyDuplicated(colnames(collision$counts)), 0L)
})
