# Regenerates fixtures/characterization.rds. Run from the package root only
# when a change of numerical results is intended (and recorded in NEWS):
#   Rscript tests/testthat/fixtures/make_characterization.R
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/characterization_helpers.R")
saveRDS(characterization_summaries(), "tests/testthat/fixtures/characterization.rds")
