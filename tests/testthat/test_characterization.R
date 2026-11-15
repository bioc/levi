library(levi)

# Characterisation test: the landscape and the regional inference must keep
# producing the same numbers while the internals are optimised. A deliberate
# change of results regenerates the fixture (see fixtures/make_characterization.R)
# and gets a "Behaviour changes" entry in NEWS.

test_that("levi() reproduces the stored landscape and regional results", {
    ref_file <- test_path("fixtures", "characterization.rds")
    skip_if_not(file.exists(ref_file), "characterization fixture missing")
    source(test_path("characterization_helpers.R"), local = TRUE)
    ref <- readRDS(ref_file)
    now <- suppressWarnings(characterization_summaries())
    expect_named(now, names(ref))
    for (nm in names(ref))
        expect_equal(now[[nm]], ref[[nm]], tolerance = 1e-8, info = nm)
})
