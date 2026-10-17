# =============================================================================
# levi -- 10. Calibration study: type I error, power and layout sensitivity
#
# Goal:     measure, by simulation, whether the permutation tests behave as
#           their documentation claims. Unit tests prove the code does what
#           it says; this script checks that the p-values are uniform under
#           the null, that the family-wise error rate stays at alpha, how
#           power grows with the effect size, and how much the landscape
#           test depends on the layout.
#
# This script is NOT run when the package is built. It takes several hours
# on four cores. It writes small CSV summaries to inst/extdata/validation/,
# which the vignette "levi_validation" loads and plots in seconds.
#
# Every test runs with the settings of the real-data analyses in the paper
# (signal_mode "logfc", logfc_k 0.7, resolution 30, smoothing 40, n_perm 999
# or exact enumeration, 100 TFCE steps), so that the calibration speaks for
# the analyses as they were run.
#
# Run with: Rscript inst/scripts/10-calibration.R [workers] [replicates]
#                   [node_replicates] [power_replicates] [node_power_replicates]
# `replicates` is the number of null replicates for the sample-label tests
# (exact 4 vs 4 enumeration, fast); `node_replicates` the number for the
# node-label test, where every replicate rebuilds 999 landscapes and is the
# cost that dominates; `power_replicates` the number per power cell for the
# sample-label and TFCE tests, and `node_power_replicates` the number per
# power cell for the node-label test.
# =============================================================================

suppressPackageStartupMessages({ library(levi); library(BiocParallel) })
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[1]) else 4L
n_rep   <- if (length(args) >= 2) as.integer(args[2]) else 1000L
n_node  <- if (length(args) >= 3) as.integer(args[3]) else 1000L
n_power <- if (length(args) >= 4) as.integer(args[4]) else 300L
n_power_node <- if (length(args) >= 5) as.integer(args[5]) else n_power
# Analysis settings shared with supplement/kegg_apriori.R and rnaseq_airway.R.
settings <- list(logfc_k = 0.7, resolutionValueInput = 30, smoothValueInput = 40)
n_perm <- 999L
n_steps <- 100L
BPPARAM <- if (workers > 1) MulticoreParam(workers) else SerialParam()
out_dir <- file.path("inst", "extdata", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
alpha <- 0.05
t_start <- Sys.time()
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

# -- Networks ----------------------------------------------------------------
# medusa: the 30-node real network shipped with levi.
# ba300: a 300-node scale-free graph (Barabasi-Albert, m = 2) laid out with
# Kamada-Kawai, closer in size and degree distribution to a STRING module.
medusa <- system.file("extdata", "medusa.dat", package = "levi")
set.seed(2026)
ba_graph <- igraph::sample_pa(300, m = 2, directed = FALSE)
igraph::V(ba_graph)$name <- paste0("g", seq_len(300))
ba <- leviFromEdges(ba_graph, layout = "kk")
networks <- list(
    medusa = list(coord = medusa, edges = NA, type = "dat",
                  genes = unique(levi:::.parseNetwork(medusa, NA, "dat")$nodes[, 1])),
    ba300  = list(coord = ba$nodes, edges = ba$edges, type = "stg",
                  genes = ba$nodes$name))

hub_module <- function(net, size = 8) {
    # The highest-degree node and its nearest neighbours: the "responding
    # module" in the power simulations.
    ni <- levi:::.networkIndex(net$coord, net$edges, net$type)
    g <- levi:::.graphFromEdges(length(ni$nodes), ni$edges)
    hub <- which.max(igraph::degree(g))
    order <- igraph::distances(g, v = hub)[1, ]
    ni$nodes[order(order)][seq_len(size)]
}

path_module <- function(net, size = 8) {
    # A peripheral module: the lowest-degree node that starts an induced
    # path of `size` nodes, grown greedily through the lowest-degree
    # unvisited neighbour. Contrasts with hub_module(): few edges, drawn at
    # the margin of the layout, little mass under the kernel.
    ni <- levi:::.networkIndex(net$coord, net$edges, net$type)
    g <- levi:::.graphFromEdges(length(ni$nodes), ni$edges)
    deg <- igraph::degree(g)
    grow <- function(start) {
        path <- start
        while (length(path) < size) {
            nb <- setdiff(igraph::neighbors(g, path[length(path)]), path)
            if (!length(nb)) return(path)
            path <- c(path, nb[which.min(deg[nb])])
        }
        path
    }
    for (start in order(deg, seq_along(deg))) {
        path <- grow(start)
        if (length(path) == size) return(ni$nodes[path])
    }
    stop("no induced path of ", size, " nodes")
}

# -- Simulation helpers ------------------------------------------------------
sim_expression <- function(genes, n_per_group = 4, effect = 0, module = NULL) {
    x <- matrix(rnorm(length(genes) * 2 * n_per_group, sd = 1), length(genes),
                dimnames = list(genes, NULL))
    groups <- rep(c("C", "T"), each = n_per_group)
    if (effect != 0 && length(module))
        x[module, groups == "T"] <- x[module, groups == "T"] + effect
    list(x = x, groups = groups)
}

regional_min_p <- function(net, lfc, n_perm = 999, smooth = 40) {
    r <- levi(expressionInput = data.frame(ID = names(lfc), logFC = lfc),
              networkCoordinatesInput = net$coord,
              networkInteractionsInput = net$edges, fileTypeInput = net$type,
              geneSymbolInput = "ID", readExpColumn = readExpColumn("logFC-logFC"),
              signal_mode = "logfc", logfc_k = settings$logfc_k,
              resolutionValueInput = settings$resolutionValueInput,
              smoothValueInput = smooth, n_perm = n_perm, .draw = FALSE)
    p <- r$regions$summary$PSpatial
    significant <- r$regions$summary$Region[r$regions$summary$Significant]
    genes <- if (length(significant)) {
        attribution <- leviRegionGenes(r, top_n = 5)
        unique(attribution$Gene[attribution$Region %in% significant])
    } else character()
    list(min_p = if (length(p)) min(p) else NA_real_, n_regions = length(p),
         sig_genes = genes)
}

replicate_min_p <- function(net, sim) {
    r <- suppressWarnings(leviReplicateInference(sim$x, sim$groups, test = "T",
        control = "C", networkCoordinatesInput = net$coord,
        networkInteractionsInput = net$edges, fileTypeInput = net$type,
        logfc_k = settings$logfc_k,
        resolutionValueInput = settings$resolutionValueInput,
        smoothValueInput = settings$smoothValueInput, n_perm = n_perm))
    p <- r$regions$summary$PSpatial
    if (length(p)) min(p) else NA_real_
}

# Node-label test on the same simulated samples: the per-gene statistic is
# the logFC (difference of group means), as in leviReplicateInference().
node_min_p <- function(net, sim) {
    lfc <- rowMeans(sim$x[, sim$groups == "T"]) - rowMeans(sim$x[, sim$groups == "C"])
    regional_min_p(net, lfc)$min_p
}

tfce_min_p <- function(net, sim) {
    r <- suppressWarnings(leviGraphTFCEInference(sim$x, sim$groups, net$coord,
        net$edges, fileTypeInput = net$type, test = "T", control = "C",
        n_perm = n_perm, n_steps = n_steps))
    min(r$statistic$PGlobal)
}

rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(...)
# Binomial summary of a vector of smallest p-values: rejection rate, its
# standard error and a 95% Wilson interval, over all simulated datasets.
# NA here means no region was detected: it is a non-rejection.
fwer_summary <- function(p) {
    n <- length(p); k <- sum(p <= alpha, na.rm = TRUE)
    z <- qnorm(0.975); est <- k / n
    centre <- (est + z^2 / (2 * n)) / (1 + z^2 / n)
    half <- z * sqrt(est * (1 - est) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    list(replicates = n, fwer = est, fwer_se = sqrt(est * (1 - est) / n),
         fwer_lo = max(0, centre - half), fwer_hi = min(1, centre + half))
}


power_summary <- function(p) {
    rejected <- !is.na(p) & p <= alpha
    power <- mean(rejected)
    list(replicates = length(p), power = power,
         power_se = sqrt(power * (1 - power) / length(p)))
}

# -- 1. Type I error under the global null ------------------------------------
# One replicate = one dataset with no effect anywhere. The family-wise error
# rate is the fraction of replicates whose smallest p-value is <= alpha; with
# regions redetected in every permutation and the maximum taken over both
# directions, it should be at most alpha.
for (nm in names(networks)) {
    net <- networks[[nm]]
    say("Type I error, node-label regional test,", nm)
    node <- bplapply(seq_len(n_node), function(i) {
        set.seed(1000 + i)
        regional_min_p(net, setNames(rnorm(length(net$genes)), net$genes))$min_p
    }, BPPARAM = BPPARAM)
    node <- unlist(node)
    add(test = "node-label regional", network = nm, design = "iid logFC",
        fwer_summary(node), alpha = alpha, no_region = mean(is.na(node)))
    say("  node-label regional", nm, "FWER", round(mean(!is.na(node) & node <= alpha), 3))
    write.csv(data.frame(test = "node-label regional", network = nm, min_p = node),
              file.path(out_dir, sprintf("null_pvalues_regional_%s.csv", nm)),
              row.names = FALSE)

    say("Type I error, sample-label replicate test,", nm)
    rep_p <- unlist(bplapply(seq_len(n_rep), function(i) {
        set.seed(2000 + i)
        replicate_min_p(net, sim_expression(net$genes))
    }, BPPARAM = BPPARAM))
    add(test = "sample-label replicate", network = nm, design = "4 vs 4, exact",
        fwer_summary(rep_p), alpha = alpha, no_region = mean(is.na(rep_p)))
    say("  sample-label replicate", nm, "FWER", round(mean(!is.na(rep_p) & rep_p <= alpha), 3))
    write.csv(data.frame(test = "sample-label replicate", network = nm, min_p = rep_p),
              file.path(out_dir, sprintf("null_pvalues_replicate_%s.csv", nm)),
              row.names = FALSE)

    say("Type I error, graph TFCE,", nm)
    tf_p <- unlist(bplapply(seq_len(n_rep), function(i) {
        set.seed(3000 + i)
        tfce_min_p(net, sim_expression(net$genes))
    }, BPPARAM = BPPARAM))
    add(test = "graph TFCE", network = nm, design = "4 vs 4, exact",
        fwer_summary(tf_p), alpha = alpha, no_region = 0)
    say("  graph TFCE", nm, "FWER", round(mean(tf_p <= alpha), 3))
    write.csv(data.frame(test = "graph TFCE", network = nm, min_p = tf_p),
              file.path(out_dir, sprintf("null_pvalues_tfce_%s.csv", nm)),
              row.names = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "type1_error.csv"),
          row.names = FALSE)

# -- 2. Power against effect size on a hub and on a peripheral module ---------
# Three tests on the same simulated samples. The node-label test sees the
# logFC of the 4 vs 4 comparison; the two sample-label tests (regional and
# TFCE) enumerate the 70 arrangements and, after the 2/n floor for
# sign-symmetric designs, both reject at the same attainable level 2/70.
# The node-label test costs 999 landscapes per replicate and runs on the
# first `n_power_node` replicates of each cell.
power_rows <- list()
test_names <- c(node = "node-label regional", replicate = "sample-label replicate",
                tfce = "graph TFCE")
for (nm in names(networks)) {
    net <- networks[[nm]]
    modules <- list(hub = hub_module(net), path = path_module(net))
    for (mod in names(modules)) {
    module <- modules[[mod]]
    for (effect in c(0.5, 1, 1.5, 2)) {
        say("Power,", nm, mod, "module, effect", effect)
        res <- bplapply(seq_len(n_power), function(i) {
            set.seed(4000 + 100 * effect + i + 10000 * (mod == "path"))
            sim <- sim_expression(net$genes, effect = effect, module = module)
            c(node = if (i <= n_power_node) node_min_p(net, sim) else NA_real_,
              replicate = replicate_min_p(net, sim), tfce = tfce_min_p(net, sim))
        }, BPPARAM = BPPARAM)
        res <- do.call(rbind, res)
        for (test in colnames(res)) {
            n_test <- if (test == "node") min(n_power_node, n_power) else n_power
            p <- res[seq_len(n_test), test]
            power_rows[[length(power_rows) + 1L]] <- data.frame(
                test = test_names[[test]],
                network = nm, module = mod, module_size = length(module),
                effect = effect, power_summary(p))
        }
    }
    }
}
write.csv(do.call(rbind, power_rows), file.path(out_dir, "power.csv"),
          row.names = FALSE)

# -- 3. Layout sensitivity ----------------------------------------------------
# Same graph, same expression, twenty force-directed layouts. How stable is
# the set of genes attributed to significant regions? Compared with TFCE,
# which never looks at coordinates and therefore gives one answer. The effect
# (a 20-gene module shifted by 3 SD, smoothing 20) is deliberately large:
# sensitivity to the layout can only be measured where the test rejects.
say("Layout sensitivity")
set.seed(5000)
module <- hub_module(networks$ba300, size = 20)
sim <- sim_expression(networks$ba300$genes, effect = 3, module = module)
lfc <- rowMeans(sim$x[, sim$groups == "T"]) - rowMeans(sim$x[, sim$groups == "C"])
layouts <- bplapply(seq_len(20), function(i) {
    set.seed(6000 + i)
    lay <- leviFromEdges(ba_graph, layout = "fr")
    net <- list(coord = lay$nodes, edges = lay$edges, type = "stg")
    out <- regional_min_p(net, lfc, smooth = 20)
    list(min_p = out$min_p, n_regions = out$n_regions, genes = out$sig_genes)
}, BPPARAM = BPPARAM)
sets <- lapply(layouts, `[[`, "genes")
jaccard <- function(a, b) if (!length(a) && !length(b)) 1 else
    length(intersect(a, b)) / length(union(a, b))
pairs <- combn(seq_along(sets), 2)
layout_rows <- data.frame(
    layout = seq_along(layouts),
    min_p = vapply(layouts, `[[`, numeric(1), "min_p"),
    n_regions = vapply(layouts, `[[`, integer(1), "n_regions"),
    n_sig_genes = lengths(sets),
    module_recovered = vapply(sets, function(s) mean(module %in% s), numeric(1)))
write.csv(layout_rows, file.path(out_dir, "layout_sensitivity.csv"),
          row.names = FALSE)
write.csv(data.frame(
    mean_pairwise_jaccard = mean(apply(pairs, 2, function(p)
        jaccard(sets[[p[1]]], sets[[p[2]]]))),
    layouts = length(sets),
    tfce_min_p = tfce_min_p(networks$ba300, sim)),
    file.path(out_dir, "layout_sensitivity_summary.csv"), row.names = FALSE)

writeLines(c(
    "# Calibration results",
    "",
    sprintf("Generated on %s by inst/scripts/10-calibration.R with %d workers,",
            format(Sys.Date()), workers),
    sprintf("%d null replicates per sample-label test, %d per node-label test,", n_rep, n_node),
    sprintf("%d (node-label: %d) replicates per power cell, in %.0f minutes.",
            n_power, n_power_node,
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
    "Settings as in the paper's analyses: signal_mode logfc, logfc_k 0.7,",
    "resolution 30, smoothing 40, n_perm 999 or exact enumeration, 100 TFCE steps.",
    "",
    "* type1_error.csv: family-wise error rate of each test under the global",
    "  null (no effect anywhere), by network, with its SE and 95% Wilson",
    "  interval over all `replicates` null datasets; no region is a non-rejection.",
    "* null_pvalues_*.csv: smallest p-value per null replicate, for QQ plots.",
    "* power.csv: probability (and SE) of at least one significant result",
    "  when an 8-gene module, the hub module or a peripheral path, is",
    "  shifted by `effect` in 4 vs 4 samples; three tests on the same data.",
    "  Datasets without regions remain in the denominator as non-rejections.",
    "* layout_sensitivity*.csv: the same data on twenty Fruchterman-Reingold",
    "  layouts of the 300-node graph; pairwise Jaccard of the genes in",
    "  significant regions, versus the layout-free TFCE result.",
    "",
    sprintf("levi %s, R %s.", packageVersion("levi"), getRversion())),
    file.path(out_dir, "README.md"))
say("Done in", round(as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
    "minutes")
