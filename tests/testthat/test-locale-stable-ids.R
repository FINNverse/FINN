library(testthat)
library(FINN)
library(data.table)

# as.factor() and sort() order through LC_COLLATE, so the same tables used to give
# one numbering under C and another under de_DE -- and the species number indexes
# the fitted parameter matrices, so par_growth[3, ] meant a different species.
test_that("resolveSiteIDs() ids do not depend on the collation locale", {
  # names chosen so that C (byte order) and de_DE (dictionary order) disagree
  sites   <- c("SB", "S_1")
  spp     <- c("Acer", "abies_x")

  tree <- rbindlist(lapply(sites, function(s) {
    dt <- toy_tree_dt()
    dt[, siteName := s]
    dt[, species_name := fifelse(species_name == "oak", spp[1], spp[2])]
    dt
  }))
  env <- CJ(siteName = sites, year = c(2000, 2010))[, `:=`(temp = 8, prec = 900)]

  ids <- function(loc) {
    old <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)
    if (is.na(suppressWarnings(Sys.setlocale("LC_COLLATE", loc)))) return(NULL)
    obs <- suppressMessages(
      makeObsData(copy(tree), plotsize = 1, aggregate_by_site = FALSE)
    )$obs_dt
    res <- suppressMessages(
      resolveSiteIDs(copy(tree), copy(env), obs, createInitCohorts = FALSE)
    )
    list(collation = list(sort(sites), sort(spp)),
         site = unique(res$siteID_dt[, .(siteName, siteID)])[order(siteName)],
         species = unique(res$species_dt)[order(species_name)])
  }

  a <- ids("C")
  b <- ids("de_DE.UTF-8")
  skip_if(is.null(a) || is.null(b), "both locales must be available")
  # guard against the test going vacuous: the two locales must really disagree
  skip_if(identical(a$collation, b$collation),
          "the two locales collate these names identically")

  expect_identical(a$site, b$site)
  expect_identical(a$species, b$species)
})
