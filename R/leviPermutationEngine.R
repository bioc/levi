# Shared machinery for the permutation tests in leviNetworkStatistics.R,
# leviSingleCellStatistics.R, leviInference.R and leviRNAseqInference.R.
#
# Every one of those tests follows the same recipe: compute a statistic on the
# observed data, recompute it on a list of null draws (permuted sample labels,
# shuffled node scores or rewired networks), keep the maximum per draw and
# compare. This file holds the recipe once. The public functions only define
# which statistic and which kind of draw they use.
#
# All randomness lives in the functions that build the draws, which run in the
# calling process. The statistic applied to each draw is deterministic, so the
# result does not depend on the BiocParallel back-end and set.seed() before a
# call reproduces it with any BPPARAM.

# -- Null draws ---------------------------------------------------------------

# Sample-label permutations, optionally within exchangeability blocks.
#
# Returns a list with $labels (list of label vectors), $exact (whether every
# distinct assignment was enumerated), $possible (their number) and $floor
# (the smallest p-value the test can report). In the exact scheme the
# observed assignment is left out of $labels; the (k + 1) / (n + 1) estimator
# adds it back.
#
# The floor is 1 / (n_perm + 1) for Monte Carlo draws. Under exact
# enumeration of a design that is balanced within every block, the
# complement arrangement (every label swapped) negates each per-gene
# statistic, and the global p-values of every test built on this engine
# (PGlobal, PSpatial) take their maximum over both directions, so the
# complement ties the observed value and the attainable floor is
# 2 / n_possible rather than 1 / n_possible. The per-direction POver and
# PUnder columns are compared with one direction's null and keep 1 / n.
.labelPermutations <- function(groups, blocks = NULL, n_perm = 999L,
                               method = c("auto", "monte_carlo", "exact"),
                               max_exact = 50000L, alpha = 0.05) {
    method <- match.arg(method)
    if (method != "exact") n_perm <- .checkNperm(n_perm)
    groups <- as.character(groups)
    if (length(unique(groups)) != 2L)
        stop("Exactly two groups are required.", call. = FALSE)
    if (is.null(blocks)) blocks <- rep("all", length(groups))
    by_block <- split(seq_along(groups), blocks)

    # Count the distinct arrangements before enumerating anything: a block
    # that holds a single condition (a donor with only controls, say) has
    # nothing to exchange and contributes one fixed arrangement. Counting
    # with choose() keeps an unblocked 58-vs-49 design (about 1e31
    # arrangements) from being expanded with combn(), which used to overflow
    # the integer range and abort before the Monte Carlo branch was reached.
    n_possible <- prod(vapply(by_block, function(i) {
        lev <- unique(groups[i])
        if (length(lev) == 1L) return(1)
        choose(length(i), sum(groups[i] == lev[2]))
    }, numeric(1)))
    exact <- method == "exact" ||
        (method == "auto" && n_possible <= max_exact)
    if (exact && n_possible > max_exact)
        stop("Exact permutation space exceeds 'max_exact'.", call. = FALSE)

    # When the floor sits above alpha no result can be declared significant,
    # which is a property of the design and not of the data, so say so up
    # front.
    n_draws <- if (exact) n_possible - 1 else as.integer(n_perm)
    symmetric <- exact && all(vapply(by_block, function(i) {
        lev <- unique(groups[i])
        length(lev) == 1L || 2L * sum(groups[i] == lev[2]) == length(i)
    }, logical(1)))
    floor_p <- if (symmetric) min(1, 2 / n_possible) else 1 / (n_draws + 1)
    if (floor_p > alpha) {
        cause <- if (symmetric) {
            sprintf(paste("Only %s distinct label arrangements are possible,",
                          "and they come in tied pairs because swapping",
                          "every label gives the same maximum statistic",
                          "over both directions"),
                    format(n_possible, big.mark = ","))
        } else if (exact) {
            sprintf("Only %s distinct label arrangements are possible",
                    format(n_possible, big.mark = ","))
        } else {
            sprintf("With n_perm = %d Monte Carlo draws", n_draws)
        }
        text <- sprintf("%s, so the smallest attainable p-value is %.3f, ",
                        cause, floor_p)
        text <- sprintf("%sabove alpha = %g. No result of this test can be ",
                        text, alpha)
        warning(text, "significant at that level.", call. = FALSE)
    }

    if (!exact) {
        labels <- replicate(as.integer(n_perm), {
            x <- groups
            for (i in by_block) x[i] <- sample(x[i])
            x
        }, simplify = FALSE)
        return(list(labels = labels, exact = FALSE, possible = n_possible,
                    floor = floor_p))
    }

    # Every way of placing the test labels inside one block.
    choices <- lapply(by_block, function(i) {
        lev <- unique(groups[i])
        if (length(lev) == 1L) return(list(groups[i]))
        test <- lev[2]
        k <- sum(groups[i] == test)
        lapply(utils::combn(seq_along(i), k, simplify = FALSE), function(pos) {
            x <- rep(lev[1], length(i))
            x[pos] <- test
            x
        })
    })
    grid <- expand.grid(lapply(choices, seq_along), KEEP.OUT.ATTRS = FALSE)
    labels <- lapply(seq_len(nrow(grid)), function(r) {
        x <- groups
        for (j in seq_along(by_block))
            x[by_block[[j]]] <- choices[[j]][[grid[r, j]]]
        x
    })
    labels <- labels[!vapply(labels, identical, logical(1), groups)]
    list(labels = labels, exact = TRUE, possible = n_possible, floor = floor_p)
}

# Number of Monte Carlo draws: a single whole number of at least 1.
.checkNperm <- function(n_perm) {
    if (!is.numeric(n_perm) || length(n_perm) != 1L || !is.finite(n_perm) ||
        n_perm < 1 || n_perm != floor(n_perm))
        stop("'n_perm' must be a single whole number >= 1.", call. = FALSE)
    as.integer(n_perm)
}

# Node-label draws: the same scores in a random order.
# Node-shuffle null that is never held in memory as a whole. Moran and
# Getis-Ord produce one statistic per node per draw, so the draws x nodes
# matrix (and the copies the p-value arithmetic makes of it) cost O(n_perm *
# n): about 1.4 GB for 10,000 nodes and 4,999 draws. Here the draws are made
# and evaluated `block` at a time and each block is folded into running
# counts by `tally(block_matrix, counts)`, so the memory is O(block * n) and
# n_perm can be raised to resolve small BH-adjusted p-values. The shuffles
# come from the same sequential calls to sample() as .nodeShuffles(), and
# the statistics use no random numbers, so the counts equal those of the
# one-piece computation.
.blockedNull <- function(x, n_perm, statistic, tally, BPPARAM = SerialParam(),
                         block = 250L) {
    n_perm <- .checkNperm(n_perm)
    counts <- NULL
    done <- 0L
    while (done < n_perm) {
        size <- min(as.integer(block), n_perm - done)
        draws <- replicate(size, sample(x), simplify = FALSE)
        counts <- tally(.permutationNull(draws, statistic, BPPARAM = BPPARAM),
                        counts)
        done <- done + size
    }
    list(counts = counts, n_perm = n_perm)
}

.nodeShuffles <- function(x, n_perm) {
    replicate(.checkNperm(n_perm), sample(x), simplify = FALSE)
}

# Degree-preserving rewiring draws: one edge index matrix per null network.
.rewiredEdgeLists <- function(edges, n_nodes, n_perm, steps) {
    g <- igraph::simplify(.graphFromEdges(n_nodes, edges))
    replicate(.checkNperm(n_perm), {
        rewired <- igraph::rewire(g,
            with = igraph::keeping_degseq(niter = steps))
        igraph::as_edgelist(rewired, names = FALSE)
    }, simplify = FALSE)
}

# -- Null distribution and p-values ------------------------------------------

# Applies `statistic` to every draw and stacks the results as rows. Runs
# through BiocParallel; the default SerialParam keeps examples and tests
# single-threaded, and MulticoreParam() or SnowParam() parallelise the loop
# without any other change.
# One-channel landscape for the permutation loops. The smoothed occupancy
# (the silhouette and the denominator of the Nadaraya-Watson average) depends
# only on the coordinates and the support weights, so it is computed once and
# every draw pays for a single blur instead of four. The result is
# bit-identical to landscape_gauss(...)$m1 for the same arguments.
.m1Landscape <- function(coord, resolution, zoom, increase, sigma, occFrac,
                         weights = numeric(0)) {
    occupancy <- landscape_occupancy(coord, resolution, zoom, increase, sigma,
                                     occFrac, weights)
    function(signal)
        landscape_m1(coord, as.numeric(signal), occupancy, resolution, zoom,
                     increase, sigma, weights)
}

# Benjamini-Hochberg across m node-wise tests needs the P-values to reach
# alpha * rank / m. With n_perm draws the smallest P-value is 1 / (n_perm + 1),
# so an isolated node can only be significant when m / (n_perm + 1) <= alpha.
# Below that resolution the family has to contain several nodes at the floor
# before anything is declared, and in practice nothing is. The warning has
# its own class so that callers can handle or silence it.
.warnResolution <- function(n_perm, m, alpha, what) {
    floor_p <- 1 / (n_perm + 1)
    if (m * floor_p <= alpha) return(invisible(FALSE))
    needed <- ceiling(m * floor_p / alpha)
    msg <- sprintf(paste0(
        "%s: with n_perm = %d the smallest attainable P-value is %.4g; across ",
        "%d nodes BH needs at least %d nodes at that floor before any is ",
        "significant at %.3g. Increase n_perm (to at least %d for a single ",
        "node to qualify)."), what, n_perm, floor_p, m, needed, alpha,
        ceiling(m / alpha) - 1L)
    warning(structure(class = c("levi_resolution_warning", "warning",
                                "condition"),
                      list(message = msg, call = NULL)))
    invisible(TRUE)
}

# 'test' and 'control' default to the order in which the group labels first
# appear, so reordering the samples silently swaps the direction of every
# statistic. The default stays for compatibility, but it now warns, and the
# resolved pair is stored in the result.
.warnGroupOrder <- function(test_missing, control_missing, test, control) {
    if (test_missing || control_missing)
        warning(structure(class = c("levi_group_order_warning", "warning",
                                    "condition"),
            list(message = sprintf(paste0(
                "'test' and 'control' were taken from the order of appearance ",
                "of the group labels (test = '%s', control = '%s'); reordering ",
                "the samples would swap the direction. Pass them explicitly."),
                test, control), call = NULL)))
    invisible(NULL)
}

# Numbers that say how much a landscape result can be trusted before any
# P-value is read: how connected the network is, and how wide the kernel is
# relative to the edges and to the picture. A kernel narrower than an edge
# smooths nothing across the network; one wider than a few edges blurs
# neighbourhoods together. None of them enters a statistic.
.networkDiagnostics <- function(node_coordinates, edge_index, sigma_coord,
                                landscape) {
    n <- nrow(node_coordinates)
    g <- .graphFromEdges(n, unique(cbind(pmin(edge_index[, 1], edge_index[, 2]),
                                         pmax(edge_index[, 1], edge_index[, 2]))))
    comp <- igraph::components(g)
    biggest <- which.max(comp$csize)
    giant <- igraph::induced_subgraph(g, which(comp$membership == biggest))
    edge_len <- if (nrow(edge_index))
        sqrt(rowSums((node_coordinates[edge_index[, 1], , drop = FALSE] -
                      node_coordinates[edge_index[, 2], , drop = FALSE])^2))
    else numeric(0)
    positive <- edge_len[edge_len > 0]
    median_len <- if (length(positive)) stats::median(positive) else NA_real_
    list(n_nodes = n, n_edges = nrow(edge_index),
         n_components = comp$no,
         largest_component_fraction = max(comp$csize) / n,
         diameter = if (igraph::vcount(giant) > 1L) igraph::diameter(giant)
                    else 0,
         median_edge_length = median_len,
         sigma_in_edge_lengths = sigma_coord / median_len,
         silhouette_fraction = mean(!is.na(landscape)))
}

# Standardised landscape field for regions. S(x) minus the mean node signal,
# divided by sd(node signals) * the layout geometry g(x) from
# landscape_null_scale(), shifted to 0.5 so that the region code can keep
# treating 0.5 as neutral. The mean and the standard deviation are those of
# the node signals of the data at hand: invariant under a node-label
# permutation, recomputed in every sample-label permutation. A constant
# signal has no spread to standardise by and gives the neutral field.
.standardiseField <- function(m1, node_signal, geometry) {
    spread <- stats::sd(node_signal)
    if (!is.finite(spread) || spread == 0) {
        m1[!is.na(m1)] <- 0.5
        return(m1)
    }
    z <- (m1 - mean(node_signal)) / (spread * geometry)
    z[!is.na(z) & !is.finite(z)] <- 0
    0.5 + z
}

# The region field of a result's metadata, as function(m1, signal).
.regionField <- function(metadata) {
    if (!identical(metadata$region_statistic, "standardised"))
        return(function(m1, signal) m1)
    function(m1, signal)
        .standardiseField(m1,
            signal[seq_along(metadata$field_nodes)][metadata$field_nodes],
            metadata$field_geometry)
}

# Significance level shared by the graph and pseudobulk tests.
.checkSigLevel <- function(sig_level) {
    if (length(sig_level) != 1L || !is.numeric(sig_level) ||
        !is.finite(sig_level) || sig_level <= 0 || sig_level >= 1)
        stop("'sig_level' must be a single number in (0, 1).", call. = FALSE)
    invisible(sig_level)
}

# Input checks shared by the permutation functions. Failing early replaces
# the silent paths found by the panel: a gene table with no gene in common
# with the network used to return zeros plus one warning per permutation,
# duplicated gene names were accepted, and an unknown 'test' level reached
# limma as an opaque error.
.checkExpressionInput <- function(genes, nodes, groups = NULL, test = NULL,
                                  control = NULL, n_perm = NULL) {
    genes <- as.character(genes)
    if (!length(genes) || is.null(genes) || anyNA(genes))
        stop("The expression matrix needs non-missing gene row names.",
             call. = FALSE)
    if (anyDuplicated(genes))
        stop("Duplicated gene names in the expression matrix (",
             sum(duplicated(genes)), " repeated, e.g. '",
             genes[duplicated(genes)][1], "'). Aggregate them first.",
             call. = FALSE)
    common <- sum(genes %in% nodes)
    if (common == 0L)
        stop("No gene of the expression matrix is a node of the network. ",
             "Check that both use the same identifiers.", call. = FALSE)
    if (common < 0.25 * length(nodes))
        warning("Only ", common, " of ", length(nodes), " network nodes ",
                "have expression data; the remaining nodes are absent from ",
                "the statistic.", call. = FALSE)
    if (!is.null(groups)) {
        if (anyNA(groups))
            stop("'groups' contains missing values.", call. = FALSE)
        if (!all(c(as.character(test), as.character(control)) %in%
                 as.character(groups)))
            stop("'test' and 'control' must both be levels of 'groups' (found: ",
                 paste(unique(as.character(groups)), collapse = ", "), ").",
                 call. = FALSE)
        if (identical(as.character(test), as.character(control)))
            stop("'test' and 'control' must differ.", call. = FALSE)
    }
    if (!is.null(n_perm) && (length(n_perm) != 1L || !is.finite(n_perm) ||
                             n_perm < 1))
        stop("'n_perm' must be a positive integer.", call. = FALSE)
    invisible(TRUE)
}

.permutationNull <- function(draws, statistic, BPPARAM = SerialParam(),
                             names = NULL) {
    rows <- bplapply(draws, statistic, BPPARAM = BPPARAM)
    null <- do.call(rbind, lapply(rows, as.numeric))
    if (is.null(null)) null <- matrix(numeric(), 0L, length(names))
    if (!is.null(names)) colnames(null) <- names
    null
}

# (k + 1) / (n + 1) upper-tail p-values of each observed value against one
# null vector of maxima. The +1 counts the observed data as a draw, so no
# finite number of permutations reports zero. A null maximum within a
# relative 1e-12 of the observed value counts as a tie, the same tolerance
# as .regionalPvalues(): under exact enumeration the complement arrangement
# reproduces the observed statistic up to floating-point rounding, and
# without the tolerance a p-value of 1/n or 2/n depended on the last bits.
.maxPvalue <- function(null_max, observed) {
    # Counting by binary search on the sorted null instead of outer() keeps
    # memory at O(n) rather than O(n x nodes): 10k nodes against 9999 draws
    # used about 2 GB.
    tolerance <- 1e-12 * pmax(1, abs(observed))
    sorted <- sort(null_max, na.last = NA)
    below <- findInterval(observed - tolerance, sorted, left.open = TRUE)
    (1 + length(sorted) - below) / (length(null_max) + 1)
}

# Largest cluster mass in each direction from a .graphRegions() summary.
.directionMaxima <- function(regions) {
    c(over  = max(c(0, regions$Mass[regions$Direction == "over"])),
      under = max(c(0, regions$Mass[regions$Direction == "under"])))
}

# -- Network helpers ----------------------------------------------------------

# Node names and a de-duplicated edge index matrix, self-loops removed.
.networkIndex <- function(networkCoordinatesInput,
                          networkInteractionsInput = NA,
                          fileTypeInput = "dat") {
    p <- .parseNetwork(networkCoordinatesInput, networkInteractionsInput,
                       fileTypeInput)
    nodes <- unique(as.character(p$nodes[, 1]))
    e <- cbind(match(as.character(p$edges[, 1]), nodes),
               match(as.character(p$edges[, 2]), nodes))
    keep <- stats::complete.cases(e) & e[, 1] != e[, 2]
    e <- e[keep, , drop = FALSE]
    # A-B and B-A are the same undirected edge; keep one of them.
    e <- unique(cbind(pmin(e[, 1], e[, 2]), pmax(e[, 1], e[, 2])))
    list(nodes = nodes, edges = e)
}

# Restricts the network to the nodes that carry a finite score and renumbers
# the edges accordingly.
.alignedGraph <- function(scores, networkCoordinatesInput,
                          networkInteractionsInput, fileTypeInput) {
    if (!is.numeric(scores) || is.null(names(scores)))
        stop("'scores' must be a numeric vector named by gene identifier.",
             call. = FALSE)
    ni <- .networkIndex(networkCoordinatesInput, networkInteractionsInput,
                        fileTypeInput)
    x <- as.numeric(scores[ni$nodes])
    keep <- is.finite(x)
    if (sum(keep) < 3L)
        stop(sum(keep), " network node(s) have a finite score; at least 3 ",
             "are needed. Check that names(scores) use the same gene ",
             "identifiers as the network.", call. = FALSE)
    map <- cumsum(keep)
    e <- ni$edges[keep[ni$edges[, 1]] & keep[ni$edges[, 2]], , drop = FALSE]
    if (nrow(e)) e <- cbind(map[e[, 1]], map[e[, 2]])
    list(nodes = ni$nodes[keep], x = x[keep], edges = e)
}

.graphFromEdges <- function(n_nodes, edges) {
    g <- igraph::make_empty_graph(n = n_nodes, directed = FALSE)
    if (nrow(edges)) g <- igraph::add_edges(g, as.vector(t(edges)))
    g
}

# Symmetric sparse adjacency matrix, so that the autocorrelation statistics
# cost O(edges) rather than O(nodes^2). `weights` gives one value per edge
# row; `self` puts ones on the diagonal (Getis-Ord Gi* includes the node
# itself). A repeated edge keeps its last weight, as an assignment would.
.adjacencyMatrix <- function(n_nodes, edges, weights = NULL, self = FALSE) {
    i <- c(edges[, 1], edges[, 2])
    j <- c(edges[, 2], edges[, 1])
    x <- if (is.null(weights)) rep(1, length(i)) else rep(weights, 2)
    if (self) {
        i <- c(i, seq_len(n_nodes))
        j <- c(j, seq_len(n_nodes))
        x <- c(x, rep(1, n_nodes))
    }
    Matrix::sparseMatrix(i = as.integer(i), j = as.integer(j), x = x,
                         dims = c(n_nodes, n_nodes), use.last.ij = TRUE)
}

# Combinatorial Laplacian D - W of a (weighted) adjacency matrix.
.laplacian <- function(W) {
    Matrix::Diagonal(x = Matrix::rowSums(W)) - W
}

# Quadratic form t(x) M x for a sparse symmetric M.
.quadraticForm <- function(x, M) {
    as.numeric(Matrix::crossprod(x, M %*% x))
}

# Eigenvectors of the `k` smallest Laplacian eigenvalues. The statistic
# needs the lowest quarter of the spectrum, and partial solvers (RSpectra,
# with or without shift-invert) were measured slower than a full symmetric
# eigen() for k = n/4, so the dense decomposition is used throughout. It is
# O(n^3): about 25 s at 2500 nodes and 3.5 min at 5000 on one core.
.lowFrequencyBasis <- function(L, k, cost_limit = 4000L) {
    n <- nrow(L)
    k <- min(k, n)
    if (n > cost_limit)
        message("Eigendecomposition of a ", n, "-node Laplacian; this takes ",
                "minutes. The spectrum test is meant for pathway-scale ",
                "networks; Moran, Getis-Ord and rewiring tests scale to ",
                "large graphs.")
    vectors <- eigen(as.matrix(L), symmetric = TRUE)$vectors
    vectors[, seq.int(n - k + 1L, n), drop = FALSE]
}

# -- Expression helpers -------------------------------------------------------

# Moderated t for the last design column: `test` versus `control`, with an
# optional block fixed effect. `trend = TRUE` is the limma-trend variant for
# logCPM-type data, where the variance depends on the mean.
.moderatedT <- function(expression, groups, levels, blocks = NULL,
                        trend = FALSE) {
    g <- factor(groups, levels = levels)
    design <- if (is.null(blocks)) stats::model.matrix(~ g) else
        stats::model.matrix(~ factor(blocks) + g)
    fit <- limma::eBayes(limma::lmFit(expression, design), trend = trend)
    stats::setNames(fit$t[, ncol(design)], rownames(expression))
}

# log2 counts per million with a pseudocount of one. With `normalize =
# "TMM"` the library sizes are scaled by edgeR's trimmed mean of M-values so
# that a few highly expressed genes do not drive the composition.
.logCPM <- function(counts, normalize = c("TMM", "none")) {
    normalize <- match.arg(normalize)
    lib <- colSums(counts)
    if (normalize == "TMM") {
        if (!requireNamespace("edgeR", quietly = TRUE))
            stop("'edgeR' is required for TMM normalisation; install it or ",
                 "use normalize = \"none\".", call. = FALSE)
        # normLibSizes() replaced calcNormFactors() in edgeR 4.0.
        tmm <- if (exists("normLibSizes", asNamespace("edgeR")))
            edgeR::normLibSizes else edgeR::calcNormFactors
        lib <- lib * tmm(counts, method = "TMM")
    }
    expr <- log2(t(t(counts) / lib * 1e6) + 1)
    rownames(expr) <- rownames(counts)
    expr
}

# Keeps the cell types whose pseudobulks cover every required column key
# (donor, or donor::condition in the paired design).
.completePseudobulkTypes <- function(pb, keys, types, paired) {
    have <- if (paired) paste(pb$donor, pb$condition, sep = "::") else pb$donor
    types[vapply(types, function(tp) all(keys %in% have[pb$cell_type == tp]),
                 logical(1))]
}

# One logCPM matrix per cell type, columns in the order of `keys`.
#
# Genes that are not expressed are dropped before normalisation and the
# limma fit: a droplet matrix carries tens of thousands of all-zero rows,
# and fitting them gives more than half of the residual variances exactly
# zero, which breaks the eBayes prior and the variance trend ("eBayes
# unreliable") for the genes that matter. A gene is kept when it has at
# least `min_count` counts in `min_samples` pseudobulks, or when it belongs
# to the network (`keep_genes`), so every node keeps a statistic.
.pseudobulkExpression <- function(pb, keys, types, paired,
                                  normalize = "TMM", keep_genes = NULL,
                                  min_count = 10L, min_samples = 2L) {
    have <- if (paired) paste(pb$donor, pb$condition, sep = "::") else pb$donor
    out <- lapply(types, function(tp) {
        at <- which(pb$cell_type == tp)
        cols <- at[match(keys, have[at])]
        counts <- pb$counts[, cols, drop = FALSE]
        expressed <- rowSums(counts >= min_count) >= min_samples
        if (!is.null(keep_genes))
            expressed <- expressed | rownames(counts) %in% keep_genes
        .logCPM(counts[expressed, , drop = FALSE], normalize)
    })
    names(out) <- types
    out
}
