library(levi)

hub_n <- system.file("extdata", "hub_network.dat",    package = "levi")
hub_e <- system.file("extdata", "hub_expression.dat", package = "levi")

# One shared levi() call — reused by all tests below (resolution=10, smooth=5)
.res <- levi(
    networkCoordinatesInput = hub_n,
    expressionInput         = hub_e,
    fileTypeInput           = "dat",
    geneSymbolInput         = "ID",
    readExpColumn           = readExpColumn("Test-Control"),
    contrastValueInput      = 50,
    resolutionValueInput    = 10,
    zoomValueInput          = 50,
    smoothValueInput        = 5
)

test_that("levi() returns a list with scores, landscape and comparison", {
    expect_type(.res, "list")
    expect_true(all(c("scores", "landscape", "comparison") %in% names(.res)))
})

test_that("scores is a data.frame with LandscapeScore in [0, 1]", {
    expect_s3_class(.res$scores, "data.frame")
    sc <- .res$scores$LandscapeScore
    expect_true(all(sc >= 0 & sc <= 1))
})

test_that("comparison field matches readExpColumn argument", {
    expect_equal(.res$comparison, "Test-Control")
})

test_that("landscape is a data.frame with numeric columns", {
    expect_s3_class(.res$landscape, "data.frame")
    expect_true(nrow(.res$landscape) > 0)
})

# --- nodes without expression value -----------------------------------------

# Until 2.0.0 this detection read fixed column indices of the edge table,
# which had drifted onto coordinate columns: the check never fired, so a total
# identifier mismatch produced a flat neutral landscape in silence.

test_that("levi(): reports how many nodes have no expression value", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    expect_message(
        levi(networkCoordinatesInput = system.file("extdata",
                 "sparse_network.dat", package = "levi"),
             expressionInput = system.file("extdata",
                 "sparse_expression.dat", package = "levi"),
             fileTypeInput  = "dat",
             geneSymbolInput = "ID",
             readExpColumn  = readExpColumn("Test-Control"),
             resolutionValueInput = 10,
             smoothValueInput     = 5),
        "10 nodes without expression value")
})

test_that("levi(): writes the unmatched node names to the log it names", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    suppressMessages(
        levi(networkCoordinatesInput = system.file("extdata",
                 "sparse_network.dat", package = "levi"),
             expressionInput = system.file("extdata",
                 "sparse_expression.dat", package = "levi"),
             fileTypeInput  = "dat",
             geneSymbolInput = "ID",
             readExpColumn  = readExpColumn("Test-Control"),
             resolutionValueInput = 10,
             smoothValueInput     = 5))

    logPath <- file.path(tempdir(), "Test-Control", "levi.log")
    expect_true(file.exists(logPath))
    expect_setequal(readLines(logPath),
                    c("G2", "G3", "G4", "G6", "G7", "G9",
                      "G10", "G11", "G13", "G14"))
})

test_that("levi(): warns when no identifier matches at all", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    # Ensembl-style ids against a network named with symbols: the landscape
    # comes back uniformly neutral, which looks exactly like "nothing changed".
    expect_warning(
        res <- levi(networkCoordinatesInput = system.file("extdata",
                        "hub_network.dat", package = "levi"),
                    expressionInput = data.frame(
                        ID      = c("ENSG001", "ENSG002"),
                        Test    = c(200, 5),
                        Control = c(10, 200)),
                    fileTypeInput  = "dat",
                    geneSymbolInput = "ID",
                    readExpColumn  = readExpColumn("Test-Control"),
                    resolutionValueInput = 10,
                    smoothValueInput     = 5),
        "None of the 9 network nodes matched")

    expect_true(all(res$scores$LandscapeScore == 0.5))
})

test_that("levi(): stays quiet when every node has a value", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    expect_no_message(
        expect_no_warning(
            levi(networkCoordinatesInput = system.file("extdata",
                     "hub_network.dat", package = "levi"),
                 expressionInput = system.file("extdata",
                     "hub_expression.dat", package = "levi"),
                 fileTypeInput  = "dat",
                 geneSymbolInput = "ID",
                 readExpColumn  = readExpColumn("Test-Control"),
                 resolutionValueInput = 10,
                 smoothValueInput     = 5)))
})

# --- early validation of the inputs ------------------------------------------

# Both checks used to surface as internal errors: indexing a missing column
# gave "undefined columns selected", and a format mismatch escaped as a dplyr
# error about seq(delimiter_edge + 1, ...). Neither named the argument to fix.

test_that("levi(): names the bad column and lists the available ones", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    err <- expect_error(
        levi(networkCoordinatesInput = system.file("extdata",
                 "hub_network.dat", package = "levi"),
             expressionInput = system.file("extdata",
                 "hub_expression.dat", package = "levi"),
             fileTypeInput  = "dat",
             geneSymbolInput = "GENE",
             readExpColumn  = readExpColumn("Test-Control"),
             resolutionValueInput = 10,
             smoothValueInput     = 5),
        "geneSymbolInput")

    expect_match(conditionMessage(err), "GENE", fixed = TRUE)
    expect_match(conditionMessage(err), "ID, Test, Control", fixed = TRUE)
})

test_that("levi(): reports a network format mismatch as such", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    hub_expr <- system.file("extdata", "hub_expression.dat", package = "levi")

    # A Medusa file declared as Pajek.
    expect_error(
        levi(networkCoordinatesInput = system.file("extdata",
                 "hub_network.dat", package = "levi"),
             expressionInput = hub_expr,
             fileTypeInput  = "net",
             geneSymbolInput = "ID",
             readExpColumn  = readExpColumn("Test-Control"),
             resolutionValueInput = 10,
             smoothValueInput     = 5),
        "no '\\*Edges' section")

    # A Pajek file declared as Medusa.
    pajek <- tempfile(fileext = ".net")
    writeLines(c("*Vertices 2", "1 \"A\" 0.1 0.2", "2 \"B\" 0.3 0.4",
                 "*Edges", "1 2 1"), pajek)

    expect_error(
        levi(networkCoordinatesInput = pajek,
             expressionInput = hub_expr,
             fileTypeInput  = "dat",
             geneSymbolInput = "ID",
             readExpColumn  = readExpColumn("Test-Control"),
             resolutionValueInput = 10,
             smoothValueInput     = 5),
        "no '\\*nodes' section")
})

test_that("levi() stores network and kernel diagnostics in the metadata", {
    d <- .res$metadata$diagnostics
    expect_named(d, c("n_nodes", "n_edges", "n_components",
                      "largest_component_fraction", "diameter",
                      "median_edge_length", "sigma_in_edge_lengths",
                      "silhouette_fraction"))
    expect_equal(d$n_nodes, 9)
    expect_equal(d$n_components, 1)
    expect_equal(d$largest_component_fraction, 1)
    expect_equal(d$diameter, 2)        # hub and eight leaves: leaf-hub-leaf
    expect_gt(d$sigma_in_edge_lengths, 0)
    expect_true(d$silhouette_fraction > 0 && d$silhouette_fraction <= 1)
})

test_that("unmeasured nodes carry no weight instead of pulling the surface to neutral", {
    net <- system.file("extdata", "medusa.dat", package = "levi")
    genes <- unique(as.character(
        levi:::.parseNetwork(net, NA, "dat")$nodes[, 1]))
    set.seed(5)
    gone <- sample(genes, 10)                    # a third of the 30 nodes
    expr <- data.frame(ID = setdiff(genes, gone), Test = 4, Control = 1)
    res <- suppressMessages(levi(expr, "dat", net, geneSymbolInput = "ID",
        readExpColumn = readExpColumn("Test-Control"),
        resolutionValueInput = 40, smoothValueInput = 20, .draw = FALSE))
    z <- res$landscape$z
    target <- 4 / (4 + 1)                        # the ratio signal of every measured node
    # Every measured point has the same value, so a weighted average over
    # measured points alone is that value everywhere. With the old neutral
    # fill the surface fell towards 0.5 around each missing gene.
    expect_equal(range(z, na.rm = TRUE), c(target, target), tolerance = 1e-9)
    w <- res$metadata$support_weights
    nodes <- res$metadata$nodes$V1
    expect_true(all(w[match(gone, nodes)] == 0))
    expect_true(all(w[seq_along(nodes)][!nodes %in% gone] == 1))
})

test_that("reciprocal edges and self-loops add no extra support in the landscape", {
    nodes <- data.frame(ID = c("A", "B", "C"), X = c(1, 2, 3), Y = c(1, 3, 1))
    messy <- data.frame(V1 = c("A", "B", "B", "B"), V2 = c("B", "A", "B", "C"))
    plain <- data.frame(V1 = c("A", "B"), V2 = c("B", "C"))
    expr <- data.frame(ID = c("A", "B", "C"), Test = c(9, 2, 1), Control = 1)
    run <- function(edges) suppressMessages(levi(expr, "stg", nodes, edges,
        geneSymbolInput = "ID", readExpColumn = readExpColumn("Test-Control"),
        resolutionValueInput = 20, smoothValueInput = 20, .draw = FALSE))
    expect_equal(nrow(run(messy)$metadata$edges), 2)
    expect_identical(run(messy)$landscape, run(plain)$landscape)
})

test_that("bandwidth_edges sets the kernel width in edge lengths", {
    net <- system.file("extdata", "medusa.dat", package = "levi")
    genes <- unique(as.character(levi:::.parseNetwork(net, NA, "dat")$nodes[, 1]))
    set.seed(3)
    expr <- data.frame(ID = genes, Test = runif(length(genes), 1, 9), Control = 1)
    run <- function(...) suppressMessages(levi(expr, "dat", net,
        geneSymbolInput = "ID", readExpColumn = readExpColumn("Test-Control"),
        resolutionValueInput = 80, .draw = FALSE, ...))
    for (b in c(0.5, 1, 2.5))
        for (zoom in c(0, 100)) {
            d <- run(bandwidth_edges = b, zoomValueInput = zoom)$metadata$diagnostics
            expect_equal(d$sigma_in_edge_lengths, b, tolerance = 0.03)
        }
    narrow <- run(bandwidth_edges = 0.5)$metadata$grid$sigma
    wide <- run(bandwidth_edges = 2.5)$metadata$grid$sigma
    expect_gt(wide, narrow)
    expect_error(run(bandwidth_edges = -1), "bandwidth_edges")
    # The default is untouched.
    expect_null(run()$metadata$bandwidth_edges)
})

test_that("region_statistic = 'standardised' uses a field of unit variance under the node-label null", {
    net <- system.file("extdata", "medusa.dat", package = "levi")
    genes <- unique(as.character(levi:::.parseNetwork(net, NA, "dat")$nodes[, 1]))
    set.seed(4)
    expr <- data.frame(ID = genes, Test = runif(length(genes), 1, 9), Control = 1)
    run <- function(...) suppressMessages(levi(expr, "dat", net,
        geneSymbolInput = "ID", readExpColumn = readExpColumn("Test-Control"),
        resolutionValueInput = 40, smoothValueInput = 30, .draw = FALSE, ...))
    default <- run()
    expect_equal(default$metadata$region_statistic, "excess")
    expect_equal(default$metadata$region_threshold, 0.1)
    std <- run(n_perm = 19, region_statistic = "standardised")
    md <- std$metadata
    expect_equal(md$region_threshold, 2)          # standard deviations
    expect_identical(std$landscape, default$landscape)   # display unchanged
    # Empirical standard deviation of the field over node-label permutations.
    coord <- rbind(md$node_coordinates, (md$node_coordinates[md$edges[, 1], ] +
                                         md$node_coordinates[md$edges[, 2], ]) / 2)
    g <- md$grid
    surface <- levi:::.m1Landscape(coord, g$resolution, g$zoom, g$increase,
                                   g$sigma, g$occupancy, md$support_weights)
    field <- levi:::.regionField(md)
    set.seed(1)
    sims <- replicate(600, {
        s <- sample(md$node_signal)
        sig <- c(s, (s[md$edges[, 1]] + s[md$edges[, 2]]) / 2)
        field(surface(sig), sig)
    })
    sdz <- apply(sims, c(1, 2), stats::sd)
    expect_gt(min(sdz, na.rm = TRUE), 0.85)
    expect_lt(max(sdz, na.rm = TRUE), 1.15)
    expect_error(run(region_statistic = "standardised", region_threshold = -1),
                 "region_threshold")
})

test_that("the analytic null scale matches permutation standard deviations", {
    set.seed(2)
    n <- 25; coord <- cbind(runif(n), runif(n))
    ed <- unique(t(apply(cbind(sample(n, 40, TRUE), sample(n, 40, TRUE)), 1, sort)))
    ed <- ed[ed[, 1] != ed[, 2], , drop = FALSE]
    allc <- rbind(coord, (coord[ed[, 1], ] + coord[ed[, 2], ]) / 2)
    w <- rep(1, n + nrow(ed)); w[2] <- 0
    w[n + which(ed[, 1] == 2 | ed[, 2] == 2)] <- 0
    occ <- landscape_occupancy(allc, 30L, -0.1, 1.2 / 29, 2.5, 0.05, w)
    scale <- levi:::landscape_null_scale(allc, matrix(as.integer(ed), ncol = 2),
                                         occ, 30L, -0.1, 1.2 / 29, 2.5, w)
    meas <- setdiff(seq_len(n), 2); s0 <- runif(n)
    sims <- replicate(1000, {
        s <- s0; s[meas] <- s0[sample(meas)]
        landscape_m1(allc, c(s, (s[ed[, 1]] + s[ed[, 2]]) / 2), occ, 30L, -0.1,
                     1.2 / 29, 2.5, w)
    })
    emp <- apply(sims, c(1, 2), stats::sd)
    ok <- !is.na(emp) & !is.na(scale)
    expect_gt(stats::cor(emp[ok], scale[ok] * stats::sd(s0[meas])), 0.99)
    expect_equal(stats::median(emp[ok] / (scale[ok] * stats::sd(s0[meas]))), 1,
                 tolerance = 0.03)
})
