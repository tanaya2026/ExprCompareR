library(testthat)

test_that("compareCancerProtein warns for invalid cancer_type input", {

  # Case 1: completely invalid cancer type
  expect_warning(
    compareCancerProtein("brain cancer"),
    "is not found in cancer_tissue_map\\$cancer"
  )

  # Case 2: slightly misspelled cancer type
  expect_warning(
    compareCancerProtein("breastcancer"),  # missing space
    "is not found in cancer_tissue_map\\$cancer"
  )

})
