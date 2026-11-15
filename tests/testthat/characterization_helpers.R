# Numerical summaries of levi() on the shipped hub and medusa networks.
# test_characterization.R compares them with fixtures/characterization.rds
# so that refactoring of the landscape and inference code cannot change a
# result unnoticed.
char_summary <- function(r) {
    out <- list(
        scores = r$scores$LandscapeScore,
        landscape = c(dim(r$landscape), colSums(r$landscape[
            vapply(r$landscape, is.numeric, logical(1))], na.rm = TRUE)))
    if (!is.null(r$regions)) {
        s <- r$regions$summary
        out$regions <- s[vapply(s, is.numeric, logical(1))]
        out$null_max_mass <- r$regions$null_max_mass
    }
    out
}

characterization_summaries <- function() {
    hub_n <- system.file("extdata", "hub_network.dat", package = "levi")
    hub_e <- system.file("extdata", "hub_expression.dat", package = "levi")
    med_n <- system.file("extdata", "medusa.dat", package = "levi")
    genes <- unique(as.character(
        levi:::.parseNetwork(med_n, NA, "dat")$nodes[, 1]))
    set.seed(1)
    med_e <- data.frame(ID = genes, Test = round(runif(length(genes), 5, 200)),
                        Control = round(runif(length(genes), 5, 200)))
    run <- function(net, expr, seed, ...) {
        set.seed(seed)
        char_summary(levi(expr, "dat", net, geneSymbolInput = "ID",
            readExpColumn = readExpColumn("Test-Control"),
            resolutionValueInput = 30, smoothValueInput = 20,
            .draw = FALSE, ...))
    }
    list(
        hub_plain  = run(hub_n, hub_e, 1, n_perm = 0),
        hub_region = run(hub_n, hub_e, 2, n_perm = 49),
        medusa_plain = run(med_n, med_e, 3, n_perm = 0),
        medusa_region = run(med_n, med_e, 4, n_perm = 49),
        medusa_degree = run(med_n, med_e, 5, n_perm = 0,
                            edge_weighting = "degree"),
        medusa_logfc = run(med_n, med_e, 6, n_perm = 49, signal_mode = "logfc",
                           logfc_k = 2))
}
