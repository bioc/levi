# Pseudobulk single-cell inference. Cells are aggregated by donor (and by
# condition in paired designs); donors, never cells, are the unit that is
# permuted, and every cell type shares the same permutation so that the joint
# maximum controls the search over cell types as well.

# Column keys shared by every cell type: donors, or donor::condition pairs in
# the paired design. Also returns the group and block vectors that go with
# those columns.
.pseudobulkDesign <- function(donors, levels, paired, condition = NULL) {
    if (paired) {
        return(list(
            keys = as.vector(outer(donors, levels, paste, sep = "::")),
            groups = rep(levels, each = length(donors)),
            blocks = rep(donors, times = length(levels))))
    }
    list(keys = donors, groups = unname(condition[donors]), blocks = NULL)
}

#' Cell-type-specific TFCE inference for pseudobulk single-cell RNA-seq
#'
#' Aggregates raw counts by donor and cell type, computes a limma-trend
#' moderated t statistic per gene and cell type on TMM-normalised log2 CPM,
#' then uses graph TFCE. Donor labels are
#' permuted jointly across all selected cell types. `PGlobal` controls the
#' search over genes, directions, regions and cell types.
#' @param counts Genes-by-cells raw count matrix.
#' @param donor,cell_type One donor and cell-type label per cell.
#' @param condition Either a named two-level vector indexed by donor (a
#'   between-donor comparison) or one two-level label per cell (a paired or
#'   repeated-condition comparison). In the latter case, labels are permuted
#'   within donor.
#' @param networkCoordinatesInput,networkInteractionsInput,fileTypeInput LEVI
#'   network inputs: node coordinates, optional interactions and the file type
#'   accepted by [levi()].
#' @param min_cells Minimum cells in each donor-cell-type pseudobulk.
#' @param normalize Library-size normalisation applied before the log2
#'   counts-per-million transform: `"TMM"` (trimmed mean of M-values via
#'   \pkg{edgeR}, the default) or `"none"` (plain counts per million).
#' @param cell_types Optional cell types to analyse. By default, types present
#' in every donor are used, ensuring joint donor permutations are valid.
#' @param n_perm Number of Monte Carlo donor-label permutations.
#' @inheritParams leviGraphTFCEInference
#' @return A list of cell-type tables, global null maxima and metadata.
#' @examples
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   set.seed(1)
#'   genes <- c("HUB", paste0("N", 1:8))
#'   # Six donors, three per condition: 20 donor-label arrangements.
#'   donor <- rep(paste0("d", 1:6), each = 4)
#'   cell_type <- rep(rep(c("A", "B"), each = 2), 6)
#'   counts <- matrix(rpois(length(genes) * length(donor), 10), length(genes),
#'     dimnames = list(genes, NULL))
#'   condition <- c(d1 = "C", d2 = "C", d3 = "C", d4 = "T", d5 = "T", d6 = "T")
#'   out <- leviSingleCellTFCEInference(counts, donor, cell_type, condition,
#'     system.file("extdata", "hub_network.dat", package = "levi"),
#'     min_cells = 2, seed = 1)
#'   out$results$A
#' }
#' @export
leviSingleCellTFCEInference <- function(counts, donor, cell_type, condition,
    networkCoordinatesInput, networkInteractionsInput = NA,
    fileTypeInput = "dat", min_cells = 10L, cell_types = NULL,
    n_perm = 999L, permutation_method = c("auto", "monte_carlo", "exact"),
    max_exact = 50000L, E = .5, H = 2, n_steps = 100L, seed = NULL,
    normalize = c("TMM", "none"), sig_level = 0.05, BPPARAM = SerialParam()) {
    .checkSigLevel(sig_level)
    normalize <- match.arg(normalize)
    if (!requireNamespace("limma", quietly = TRUE))
        stop("'limma' is required.", call. = FALSE)
    paired <- length(condition) == ncol(counts)
    if (!paired && is.null(names(condition)))
        stop("A donor-level 'condition' must be named by donor.",
             call. = FALSE)
    if (length(unique(condition)) != 2L)
        stop("'condition' must have two levels.", call. = FALSE)
    if (!is.null(seed)) withr::local_seed(seed)

    pb <- if (paired) {
        leviPseudobulk(counts, donor, cell_type, condition = condition,
                       min_cells = min_cells)
    } else {
        leviPseudobulk(counts, donor, cell_type, min_cells = min_cells)
    }
    donors <- if (paired) sort(unique(as.character(donor))) else
        names(condition)
    levels <- unique(as.character(condition))
    design <- .pseudobulkDesign(donors, levels, paired, condition)

    types <- unique(pb$cell_type)
    if (!is.null(cell_types)) types <- intersect(types, cell_types)
    types <- .completePseudobulkTypes(pb, design$keys, types, paired)
    if (!length(types))
        stop("No requested cell type has pseudobulk for every donor.",
             call. = FALSE)
    ni <- .networkIndex(networkCoordinatesInput, networkInteractionsInput,
                        fileTypeInput)
    .checkExpressionInput(rownames(counts), ni$nodes, n_perm = n_perm)
    expression <- .pseudobulkExpression(pb, design$keys, types, paired,
                                        normalize,
                                        keep_genes = ni$nodes)

    statistic <- function(expr, g)
        .moderatedT(expr, g, levels, design$blocks, trend = TRUE)
    observed <- lapply(expression, function(expr) {
        t_stat <- statistic(expr, design$groups)
        list(t = t_stat,
             tfce = .tfce(t_stat, ni$nodes, ni$edges, E, H, n_steps))
    })

    # One draw yields the over/under TFCE maxima of every cell type; the
    # global null is the largest of them.
    perms <- .labelPermutations(design$groups, design$blocks, n_perm,
                                permutation_method, max_exact, alpha = sig_level)
    null <- .permutationNull(perms$labels, function(g) {
        unlist(lapply(expression, function(expr)
            .tfceMaxima(statistic(expr, g), ni$nodes, ni$edges, E, H,
                        n_steps)))
    }, BPPARAM = BPPARAM,
       names = paste(rep(types, each = 2), c("over", "under"), sep = "."))
    null_global <- apply(null, 1L, max)

    tables <- lapply(seq_along(types), function(k) {
        o <- observed[[k]]
        over_col <- 2L * k - 1L
        p_global <- .maxPvalue(null_global, pmax(o$tfce[, 1], o$tfce[, 2]))
        data.frame(
            Gene = ni$nodes, T = as.numeric(o$t[ni$nodes]),
            TFCEOver = o$tfce[, 1], TFCEUnder = o$tfce[, 2],
            POver = .maxPvalue(null[, over_col], o$tfce[, 1]),
            PUnder = .maxPvalue(null[, over_col + 1L], o$tfce[, 2]),
            PGlobal = p_global, GlobalSignificant = p_global <= sig_level)
    })
    names(tables) <- types

    list(results = tables, null_global = null_global, exact = perms$exact,
         possible_permutations = perms$possible, p_floor = perms$floor,
         donors = donors,
         cell_types = types,
         method = if (paired) {
             paste("paired donor-by-condition pseudobulk limma t; graph TFCE;",
                   "joint within-donor label maximum")
         } else {
             paste("pseudobulk donor-level limma t; graph TFCE;",
                   "joint donor-label maximum")
         })
}

#' Run Moran, Laplacian/Fourier and rewiring tests per cell type
#'
#' Takes the moderated t-statistics stored in a [leviSingleCellTFCEInference()]
#' result and applies [leviGraphMoran()], [leviGraphSpectrum()] and
#' [leviGraphRewiringInference()] to every cell type.
#' @param single_cell_tfce Result of [leviSingleCellTFCEInference()]; the `T`
#'   column of every cell-type table supplies the gene scores.
#' @param networkCoordinatesInput,networkInteractionsInput,fileTypeInput LEVI
#'   network inputs: node coordinates, optional interactions and the file type
#'   accepted by [levi()].
#' @param threshold Absolute t-statistic threshold for cluster formation in the
#'   rewiring test.
#' @param n_perm Number of node-label permutations or rewired networks used by
#'   each topology test.
#' @param seed Optional RNG seed.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling how
#'   the permutations are evaluated. The default runs them serially.
#' @return A named list of topology-test results, one entry for each cell type.
#' @examples
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   set.seed(1)
#'   net <- system.file("extdata", "hub_network.dat", package = "levi")
#'   genes <- c("HUB", paste0("N", 1:8))
#'   donor <- rep(paste0("d", 1:6), each = 4)
#'   cell_type <- rep(rep(c("A", "B"), each = 2), 6)
#'   counts <- matrix(rpois(length(genes) * length(donor), 10), length(genes),
#'     dimnames = list(genes, NULL))
#'   condition <- c(d1 = "C", d2 = "C", d3 = "C", d4 = "T", d5 = "T", d6 = "T")
#'   tfce <- leviSingleCellTFCEInference(counts, donor, cell_type, condition,
#'     net, min_cells = 2, seed = 1)
#'   topology <- leviSingleCellTopologyInference(tfce, net, threshold = 1,
#'     n_perm = 19, seed = 1)
#'   topology$A$moran$global
#' }
#' @export
leviSingleCellTopologyInference <- function(single_cell_tfce,
    networkCoordinatesInput, networkInteractionsInput = NA,
    fileTypeInput = "dat", threshold = 2, n_perm = 999L, seed = NULL,
    BPPARAM = SerialParam()) {
    lapply(single_cell_tfce$results, function(tab) {
        s <- stats::setNames(tab$T, tab$Gene)
        list(
            moran = leviGraphMoran(s, networkCoordinatesInput,
                networkInteractionsInput, fileTypeInput, n_perm, seed,
                BPPARAM = BPPARAM),
            spectrum = leviGraphSpectrum(s, networkCoordinatesInput,
                networkInteractionsInput, fileTypeInput, n_perm, seed,
                BPPARAM = BPPARAM),
            rewiring = leviGraphRewiringInference(s, networkCoordinatesInput,
                networkInteractionsInput, fileTypeInput, threshold, n_perm,
                seed = seed, BPPARAM = BPPARAM))
    })
}

#' Joint regional inference for paired single-cell pseudobulks
#'
#' Computes graph-connected cluster mass and explicit two-dimensional LEVI areas
#' for every selected cell type. Treatment labels are swapped within donors and
#' the null stores the maximum region mass across cell types, directions and
#' regions, controlling the complete search family.
#'
#' @param counts Genes-by-cells raw count matrix with gene row names.
#' @param donor,cell_type One donor and cell-type label per cell.
#' @param condition Two-level condition label per cell; every donor must
#'   contribute cells to both levels within each analysed cell type.
#' @param networkCoordinatesInput,networkInteractionsInput,fileTypeInput LEVI
#'   network inputs: node coordinates, optional interactions and the file type
#'   accepted by [levi()].
#' @param min_cells Minimum cells in each donor-condition-cell-type pseudobulk.
#' @param normalize Library-size normalisation applied before the log2
#'   counts-per-million transform: `"TMM"` (trimmed mean of M-values via
#'   \pkg{edgeR}, the default) or `"none"` (plain counts per million).
#' @param cell_types Optional cell types to analyse. By default, types with a
#'   complete donor-by-condition pseudobulk set are used.
#' @param threshold Absolute moderated-t threshold for graph cluster formation.
#' @param n_perm Number of Monte Carlo within-donor label permutations.
#' @inheritParams leviGraphTFCEInference
#' @param ... Landscape arguments passed to [levi()] (for example
#'   `resolutionValueInput` or `region_threshold`). Do not supply expression,
#'   network, signal mode, permutation or inference arguments. `logfc_k` is
#'   applied to the observed and to every permuted landscape.
#' @return A list of graph-cluster and two-dimensional landscape regional tests.
#' @examples
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   set.seed(1)
#'   genes <- c("HUB", paste0("N", 1:8))
#'   # Five paired donors: 32 within-donor swaps; 19 Monte Carlo draws here
#'   # keep the example fast.
#'   donor <- rep(paste0("d", 1:5), each = 4)
#'   condition <- rep(rep(c("ctrl", "stim"), each = 2), 5)
#'   cell_type <- rep(c("A", "B"), times = 10)
#'   counts <- matrix(rpois(length(genes) * length(donor), 20), length(genes),
#'     dimnames = list(genes, NULL))
#'   out <- leviSingleCellRegionalInference(counts, donor, cell_type, condition,
#'     system.file("extdata", "hub_network.dat", package = "levi"),
#'     min_cells = 1, threshold = 1, n_perm = 19,
#'     permutation_method = "monte_carlo", seed = 1,
#'     resolutionValueInput = 1, region_threshold = .02)
#'   out$cluster_mass$A$summary
#' }
#' @export
leviSingleCellRegionalInference <- function(counts, donor, cell_type,
    condition, networkCoordinatesInput, networkInteractionsInput = NA,
    fileTypeInput = "dat", min_cells = 10L, cell_types = NULL, threshold = 2,
    n_perm = 999L, permutation_method = c("auto", "monte_carlo", "exact"),
    max_exact = 50000L, seed = NULL, normalize = c("TMM", "none"),
    sig_level = 0.05, BPPARAM = SerialParam(), ...) {
    .checkSigLevel(sig_level)
    normalize <- match.arg(normalize)
    if (!requireNamespace("limma", quietly = TRUE))
        stop("'limma' is required.", call. = FALSE)
    if (length(condition) != ncol(counts) || length(unique(condition)) != 2L)
        stop("'condition' must be a two-level, per-cell vector.",
             call. = FALSE)
    if (!is.null(seed)) withr::local_seed(seed)

    pb <- leviPseudobulk(counts, donor, cell_type, condition = condition,
                         min_cells = min_cells)
    donors <- sort(unique(as.character(donor)))
    levels <- unique(as.character(condition))
    design <- .pseudobulkDesign(donors, levels, paired = TRUE)
    types <- unique(pb$cell_type)
    if (!is.null(cell_types)) types <- intersect(types, cell_types)
    types <- .completePseudobulkTypes(pb, design$keys, types, paired = TRUE)
    if (!length(types))
        stop("No cell type contains every donor-condition pseudobulk.",
             call. = FALSE)
    ni <- .networkIndex(networkCoordinatesInput, networkInteractionsInput,
                        fileTypeInput)
    .checkExpressionInput(rownames(counts), ni$nodes, n_perm = n_perm)
    expression <- .pseudobulkExpression(pb, design$keys, types, paired = TRUE,
                                        normalize,
                                        keep_genes = ni$nodes)

    statistic <- function(expr, g)
        .moderatedT(expr, g, levels, design$blocks, trend = TRUE)
    perms <- .labelPermutations(design$groups, design$blocks, n_perm,
                                permutation_method, max_exact, alpha = sig_level)

    # 1. Graph-connected cluster mass, one joint null over cell types.
    observed_clusters <- lapply(expression, function(expr)
        list(summary = .graphRegions(statistic(expr, design$groups),
                                     ni$nodes, ni$edges, threshold)))
    null_cluster <- .permutationNull(perms$labels, function(g) {
        per_type <- vapply(expression, function(expr)
            .directionMaxima(.graphRegions(statistic(expr, g), ni$nodes,
                                           ni$edges, threshold)),
            numeric(2))
        apply(per_type, 1L, max)
    }, BPPARAM = BPPARAM, names = c("over", "under"))
    clusters <- lapply(observed_clusters, function(z) {
        z <- .regionalPvalues(z, null_cluster)
        z$summary$Significant <- z$summary$PSpatial <= sig_level
        z
    })

    # 2. Two-dimensional landscape regions on the mean logFC per cell type.
    mean_lfc <- function(expr, g) {
        rowMeans(expr[, g == levels[2], drop = FALSE]) -
            rowMeans(expr[, g == levels[1], drop = FALSE])
    }
    parsed <- .parseNetwork(networkCoordinatesInput, networkInteractionsInput,
                            fileTypeInput)
    run_landscape <- function(expr) {
        lfc <- mean_lfc(expr, design$groups)
        do.call(levi, c(list(
            expressionInput = data.frame(.Gene = rownames(expr), logFC = lfc),
            geneSymbolInput = ".Gene",
            readExpColumn = readExpColumn("logFC-logFC"),
            signal_mode = "logfc", n_perm = 0L, inference_unit = "region",
            plot3d = FALSE,
            networkCoordinatesInput = networkCoordinatesInput,
            networkInteractionsInput = networkInteractionsInput,
            fileTypeInput = fileTypeInput,
            .parsed_network = parsed, .draw = FALSE), list(...)))
    }
    landscapes <- lapply(expression, run_landscape)

    # The grid, coordinates and edges are the same for every cell type.
    first <- landscapes[[1]]$metadata
    grid <- first$grid
    edges <- first$edges
    node_names <- as.character(first$nodes[, 1])
    coord <- rbind(first$node_coordinates,
        (first$node_coordinates[edges[, 1], , drop = FALSE] +
         first$node_coordinates[edges[, 2], , drop = FALSE]) / 2)
    # The null signal must use the sigmoid steepness of the observed landscape.
    k <- list(...)[["logfc_k"]] %||% 1
    landscape_of <- .m1Landscape(coord, grid$resolution, grid$zoom,
        grid$increase, grid$sigma, grid$occupancy,
        first$support_weights %||% numeric(0))
    null_area <- .permutationNull(perms$labels, function(g) {
        per_type <- vapply(types, function(tp) {
            node_lfc <- unname(mean_lfc(expression[[tp]], g)[node_names])
            signal <- .networkSignals(cbind(node_lfc, node_lfc), edges, TRUE,
                                      "logfc", k)$signal
            z <- .regionField(landscapes[[tp]]$metadata)(landscape_of(signal),
                                                         signal)
            .maximumRegionMass(z, grid$increase,
                landscapes[[tp]]$metadata$region_threshold,
                landscapes[[tp]]$metadata$region_min_cells)
        }, numeric(2))
        apply(per_type, 1L, max)
    }, BPPARAM = BPPARAM, names = c("over", "under"))
    landscapes <- lapply(landscapes, function(z) {
        z$regions <- .regionalPvalues(z$regions, null_area)
        z$regions$summary$Significant <- z$regions$summary$PSpatial <= sig_level
        z$regions$inference <- list(method = "joint donor-label maximum",
                                    scope = "cell types and directions")
        .annotateReplicateRegions(z)
    })

    list(cluster_mass = clusters, landscapes = landscapes,
         null_cluster = null_cluster, null_area = null_area,
         exact = perms$exact, possible_permutations = perms$possible,
         p_floor = perms$floor, donors = donors, cell_types = types,
         method = paste("paired pseudobulk; joint within-donor maximum",
                        "across cell types"))
}

#' Cell-type-specific condition interaction with joint graph TFCE inference
#'
#' Tests whether the condition effect in each non-reference cell type differs
#' from the reference type. For every donor the interaction contrast
#' \eqn{(y_{stim} - y_{ctrl})_{type} - (y_{stim} - y_{ctrl})_{reference}} is
#' formed on TMM log2 CPM pseudobulks and tested with a limma-trend
#' one-sample moderated t. A condition effect shared by all cell types
#' cancels in the contrast, so swapping the condition labels inside a donor
#' (jointly over all cell types) flips the sign of that donor's contrast and
#' is a valid null for "no interaction" even when the main effect is large.
#' Cell types are aligned on the genes they share. Every donor must have
#' both conditions in every analysed cell type.
#'
#' @inheritParams leviSingleCellRegionalInference
#' @param reference_cell_type Cell type against which the condition effect of
#'   every other type is contrasted. Defaults to the first complete type.
#' @param cell_types Optional cell types to analyse; at least two complete
#'   types are required.
#' @inheritParams leviGraphTFCEInference
#' @return A list of interaction TFCE tables and the joint donor-label null
#'   distribution.
#' @examples
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   set.seed(1)
#'   genes <- c("HUB", paste0("N", 1:8))
#'   donor <- rep(paste0("d", 1:5), each = 4)
#'   condition <- rep(rep(c("ctrl", "stim"), each = 2), 5)
#'   cell_type <- rep(c("A", "B"), times = 10)
#'   counts <- matrix(rpois(length(genes) * length(donor), 20), length(genes),
#'     dimnames = list(genes, NULL))
#'   out <- leviSingleCellInteractionTFCE(counts, donor, cell_type, condition,
#'     system.file("extdata", "hub_network.dat", package = "levi"),
#'     reference_cell_type = "A", min_cells = 1, n_perm = 19,
#'     permutation_method = "monte_carlo", seed = 1)
#'   out$results$B
#' }
#' @export
leviSingleCellInteractionTFCE <- function(counts, donor, cell_type, condition,
    networkCoordinatesInput, networkInteractionsInput = NA,
    fileTypeInput = "dat", reference_cell_type = NULL, min_cells = 10L,
    cell_types = NULL, n_perm = 999L,
    permutation_method = c("auto", "monte_carlo", "exact"),
    max_exact = 50000L, E = .5, H = 2, n_steps = 100L, seed = NULL,
    normalize = c("TMM", "none"), sig_level = 0.05, BPPARAM = SerialParam()) {
    .checkSigLevel(sig_level)
    normalize <- match.arg(normalize)
    if (!requireNamespace("limma", quietly = TRUE))
        stop("'limma' is required.", call. = FALSE)
    if (length(condition) != ncol(counts) || length(unique(condition)) != 2L)
        stop("'condition' must be a two-level label per cell.", call. = FALSE)
    if (!is.null(seed)) withr::local_seed(seed)

    pb <- leviPseudobulk(counts, donor, cell_type, condition = condition,
                         min_cells = min_cells)
    donors <- sort(unique(as.character(donor)))
    levels <- unique(as.character(condition))
    design <- .pseudobulkDesign(donors, levels, paired = TRUE)
    types <- unique(pb$cell_type)
    if (!is.null(cell_types)) types <- intersect(types, cell_types)
    types <- .completePseudobulkTypes(pb, design$keys, types, paired = TRUE)
    if (length(types) < 2L)
        stop("At least two complete cell types are required.", call. = FALSE)
    if (is.null(reference_cell_type)) reference_cell_type <- types[1]
    if (!reference_cell_type %in% types)
        stop("Reference cell type is unavailable.", call. = FALSE)
    others <- setdiff(types, reference_cell_type)
    ni <- .networkIndex(networkCoordinatesInput, networkInteractionsInput,
                        fileTypeInput)
    .checkExpressionInput(rownames(counts), ni$nodes, n_perm = n_perm)
    expression <- .pseudobulkExpression(pb, design$keys, types, paired = TRUE,
                                        normalize,
                                        keep_genes = ni$nodes)

    # Per-donor interaction contrasts. For donor j and cell type tp the
    # condition effect is diff_tp[, j] = logCPM(stim) - logCPM(ctrl); the
    # interaction is d_tp[, j] = diff_tp[, j] - diff_ref[, j]. A common
    # condition effect cancels in d, so under the null of no interaction d
    # is symmetric about zero and swapping the condition labels inside a
    # donor (jointly over all cell types) only flips the sign of that donor's
    # column. The statistic is a one-sample moderated t on d. The earlier
    # version refitted `~ donor + type * condition` on the swapped labels,
    # which is invariant only under the global null (no condition effect at
    # all): with a common main effect and no interaction it rejected in 85%
    # of null data sets (expert panel, 2026-11). The cell types are aligned
    # on their common genes, because .pseudobulkExpression() filters each
    # type separately and stacking them by position mixed genes.
    common <- Reduce(intersect, lapply(expression, rownames))
    if (!length(intersect(ni$nodes, common)))
        stop("No network node is shared by the selected cell types.",
             call. = FALSE)
    expression <- lapply(expression, function(e) e[common, , drop = FALSE])
    first <- design$groups == levels[1]
    second <- design$groups == levels[2]
    donor_order <- design$blocks[second]
    if (!identical(design$blocks[first], donor_order))
        stop("Internal error: paired pseudobulk columns are misaligned.",
             call. = FALSE)
    effect <- lapply(expression, function(e)
        e[, second, drop = FALSE] - e[, first, drop = FALSE])
    contrasts <- lapply(others, function(tp)
        effect[[tp]] - effect[[reference_cell_type]])
    names(contrasts) <- others
    abundance <- rowMeans(do.call(cbind, expression))

    # Sign of each donor's contrast under a label vector: +1 when the donor
    # keeps its observed labels, -1 when they are swapped.
    donorSigns <- function(g) {
        ifelse(g[second] == levels[2], 1, -1)
    }
    statistic <- function(g) {
        s <- donorSigns(g)
        out <- lapply(contrasts, function(d) {
            fit <- limma::lmFit(t(t(d) * s))
            fit$Amean <- abundance
            fit <- limma::eBayes(fit, trend = TRUE)
            stats::setNames(fit$t[, 1], common)
        })
        names(out) <- others
        out
    }

    observed_t <- statistic(design$groups)
    observed <- lapply(observed_t, function(s)
        .tfce(s, ni$nodes, ni$edges, E, H, n_steps))

    perms <- .labelPermutations(design$groups, design$blocks, n_perm,
                                permutation_method, max_exact, alpha = sig_level)
    null <- .permutationNull(perms$labels, function(g) {
        max(vapply(statistic(g), function(s)
            max(.tfceMaxima(s, ni$nodes, ni$edges, E, H, n_steps)),
            numeric(1)))
    }, BPPARAM = BPPARAM, names = "global")
    null_global <- null[, "global"]

    tables <- lapply(others, function(tp) {
        o <- observed[[tp]]
        p_global <- .maxPvalue(null_global, pmax(o[, 1], o[, 2]))
        data.frame(
            Gene = ni$nodes, T = as.numeric(observed_t[[tp]][ni$nodes]),
            TFCEOver = o[, 1], TFCEUnder = o[, 2], PGlobal = p_global,
            GlobalSignificant = p_global <= sig_level,
            Interaction = paste0(tp, " versus ", reference_cell_type))
    })
    names(tables) <- others

    list(results = tables, null_global = null_global, exact = perms$exact,
         possible_permutations = perms$possible, p_floor = perms$floor,
         reference_cell_type = reference_cell_type,
         method = paste("pseudobulk per-donor interaction contrasts,",
                        "moderated one-sample t; joint within-donor",
                        "sign-flip TFCE maximum"))
}
