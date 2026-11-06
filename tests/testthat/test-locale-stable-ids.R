test_that("site, patch and species ids do not depend on the collation locale", {
  skip_if_no_torch()
  # as.factor() and sort() order through LC_COLLATE, so the same file used to give one
  # numbering under C and another under de_DE -- and the species number indexes the
  # fitted parameter matrices, so par_growth[3, ] meant a different species.
  f <- system.file("extdata", "example_tree_dt.csv", package = "FINN")
  skip_if(f == "", "example_tree_dt.csv not installed")
  tree <- data.table::fread(f)
  ids <- function(loc) {
    old <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)
    if (is.na(suppressWarnings(Sys.setlocale("LC_COLLATE", loc)))) return(NULL)
    s <- sort(unique(tree$siteName), method = "radix")
    list(site = match(unique(tree$siteName), s),
         species = match(unique(tree$species_name),
                         sort(unique(tree$species_name), method = "radix")))
  }
  a <- ids("C")
  b <- ids("de_DE.UTF-8")
  skip_if(is.null(a) || is.null(b), "both locales must be available")
  expect_identical(a$site, b$site)
  expect_identical(a$species, b$species)
})
