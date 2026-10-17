# Evaluate summary helpers only; never run the simulation study during tests.
calibration_summaries <- function() {
    script <- system.file("scripts", "10-calibration.R", package = "levi")
    env <- new.env(parent = baseenv())
    env$alpha <- .05
    env$qnorm <- stats::qnorm
    for (expr in parse(script)) {
        if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
            as.character(expr[[2]]) %in% c("fwer_summary", "power_summary"))
            eval(expr, env)
    }
    env
}

test_that("calibration counts datasets without regions as non-rejections", {
    env <- calibration_summaries()
    p <- c(.01, NA_real_, .5, NA_real_)
    fwer <- env$fwer_summary(p)
    power <- env$power_summary(p)
    expect_equal(fwer$replicates, 4L)
    expect_equal(fwer$fwer, .25)
    expect_equal(fwer$fwer_se, sqrt(.25 * .75 / 4))
    expect_equal(power$replicates, 4L)
    expect_equal(power$power, .25)
    expect_equal(power$power_se, fwer$fwer_se)
    expect_equal(env$fwer_summary(rep(NA_real_, 3))$fwer, 0)
    expect_equal(env$power_summary(rep(NA_real_, 3))$power, 0)
})
