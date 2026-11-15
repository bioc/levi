library(levi)

# The C++ maximum regional mass must agree with the R reference.

test_that("C++ maximum region mass matches the R reference", {
    set.seed(11)
    for (rep in 1:80) {
        nr <- sample(c(3, 12, 40), 1); nc <- sample(c(3, 12, 40), 1)
        z <- matrix(0.5 + rnorm(nr * nc, sd = 0.2), nr, nc)
        z[sample(length(z), max(1, length(z) %/% 10))] <- NA
        if (rep %% 5 == 0) z[] <- 0.5
        thr <- sample(c(0, 0.05, 0.1, 0.2), 1)
        mc <- sample(c(1L, 3L, 10L), 1)
        inc <- runif(1, 0.01, 1)
        expect_equal(levi:::.maximumRegionMass(z, inc, thr, mc),
                     levi:::.maximumRegionMassR(z, inc, thr, mc),
                     tolerance = 1e-12,
                     info = sprintf("rep %d", rep))
    }
})

test_that("the blocked node-shuffle null does not depend on the block size", {
    x <- c(3, 1, 4, 1, 5, 9, 2, 6)
    stat <- function(o) cumsum(x[o])
    obs <- cumsum(x)
    tally <- function(b, counts) {
        hit <- colSums(b >= rep(obs, each = nrow(b)))
        if (is.null(counts)) hit else counts + hit
    }
    run <- function(block) { set.seed(5)
        levi:::.blockedNull(seq_along(x), 103, stat, tally, block = block) }
    expect_identical(run(7L)$counts, run(1000L)$counts)
    expect_identical(run(1L)$counts, run(1000L)$counts)
    expect_equal(run(7L)$n_perm, 103L)
})
