library(testthat)

test_that("detect_outliers warns for invalid tissue input", {

  # Case 1: completely invalid tissue
  expect_warning(
    detect_outliers("brainstem"),
    "is not found in tissue_map\\$protein_tissue"
  )

  # Case 2: slightly misspelled tissue
  expect_warning(
    detect_outliers("lvg"),  # typo instead of "liver"
    "is not found in tissue_map\\$protein_tissue"
  )

})
