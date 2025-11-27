library(testthat)

# This file includes tests for function detect_outliers

#-------------------------------------------------------------------------------
# Unit Tests

# Testing detect_outliers for valid inputs

test_that("detect_outliers warns for invalid tissue input", {

  # Case 1: completely invalid tissue
  expect_error(
    detect_outliers("brainstem"),
    "is not found in tissue_map\\$protein_tissue"
  )

  # Case 2: slightly misspelled tissue
  expect_error(
    detect_outliers("lvg"),  # typo instead of "liver"
    "is not found in tissue_map\\$protein_tissue"
  )

})

#-------------------------------------------------------------------------------
# Integration Tests

# Test if detect_outliers returns a vector and plot; the vector must be a character vector
# Also test that for a tissue like 'endometrium' there are no outliers
test_that("detect_outliers returns empty outlier_vector for endometrium", {

  output <- detect_outliers("endometrium")

  # Check the structure
  expect_true(is.list(output))
  expect_true("outlier_vector" %in% names(output))
  expect_true("outlier_plot" %in% names(output))

  # outlier_vector should be an empty character vector
  expect_type(output$outlier_vector, "character")
  expect_length(output$outlier_vector, 0)

  # outlier_plot should be a plot (ggplot)
  expect_true(is.object(output$outlier_plot))
})

#-------------------------------------------------------------------------------

# Integration Tests

# Test if detect_outliers returns a vector and plot; the vector must be a character vector
# Also test that detect_outliers returns the correct vector with the correct outliers.


test_that("detect_outliers returns correct vector and plot for lung", {

  output <- detect_outliers("lung")

  # Check structure
  expect_true(is.list(output))
  expect_true(all(c("outlier_vector", "outlier_plot") %in% names(output)))

  # Expected outlier genes
  expected_genes <- c(
    "CD74", "VIM", "ENO1", "ACTB", "HSP90AB1", "HSPB1", "ANXA10", "SPARC", "PEX5L", "FN1",
    "AGMAT", "NT5C1A", "CTSD", "SFTPA1", "SRGN", "LIN28A", "FCRLA", "ZNF593", "RPS27A", "MED12L",
    "GC", "SPINK7", "SPHKAP", "ADCY8", "EEF1A1", "GRHL3", "RPL8", "S100A11", "S100A9", "BSN",
    "ZMAT4", "SFTPC", "SFTPB", "PURG", "SFTPA2", "GJB4", "HLA-DRB1", "S100A4", "PSAP", "S100A6",
    "HLA-DRA", "AGER", "HLA-E", "FAM221B", "XKR9", "HLA-B", "AQP1", "HBB", "EEF1G"
  )

  # Check vector type and contents
  expect_type(output$outlier_vector, "character")
  expect_equal(output$outlier_vector, expected_genes)

  # Check plot type
  expect_true(is.object(output$outlier_plot))
})

#-------------------------------------------------------------------------------

# [END]

