library(levi)

# The C++ TFCE (src/tfce.cpp, union-find) must agree with the plain-R
# reference .tfceR() on graphs with isolated nodes, tied statistics, missing
# values and degenerate thresholds.

random_graph <- function(n, m) {
    e <- if (m > 0) cbind(sample(n, m, TRUE), sample(n, m, TRUE)) else
        matrix(integer(), 0L, 2L)
    e <- e[e[, 1] != e[, 2], , drop = FALSE]
    unique(cbind(pmin(e[, 1], e[, 2]), pmax(e[, 1], e[, 2])))
}

test_that("C++ TFCE matches the R reference", {
    set.seed(42)
    for (rep in 1:60) {
        n <- sample(c(5, 30, 120), 1)
        e <- random_graph(n, sample(c(0, n, 3 * n), 1))
        nodes <- paste0("g", seq_len(n))
        x <- stats::setNames(
            if (rep %% 3 == 0) round(rnorm(n), 1) else rnorm(n), nodes)  # ties
        if (rep %% 7 == 0) x[sample(n, 2)] <- NA
        if (rep %% 11 == 0) x[] <- 0
        steps <- sample(c(1L, 2L, 5L, 100L), 1)
        fast <- levi:::.tfce(x, nodes, e, 0.5, 2, steps)
        slow <- suppressWarnings(levi:::.tfceR(x, nodes, e, 0.5, 2, steps))
        expect_equal(fast, slow, tolerance = 1e-12,
                     info = sprintf("rep %d (n = %d, steps = %d)", rep, n, steps))
    }
})

test_that("C++ TFCE honours E and H and returns zeros without signal", {
    nodes <- paste0("g", 1:6)
    e <- cbind(1:5, 2:6)
    x <- stats::setNames(c(2, 3, 1, -2, -3, -1), nodes)
    for (EH in list(c(0.5, 2), c(1, 2), c(0, 1), c(2, 0.5)))
        expect_equal(levi:::.tfce(x, nodes, e, EH[1], EH[2], 50L),
                     levi:::.tfceR(x, nodes, e, EH[1], EH[2], 50L),
                     tolerance = 1e-12)
    flat <- stats::setNames(rep(0, 6), nodes)
    expect_true(all(levi:::.tfce(flat, nodes, e) == 0))
})
