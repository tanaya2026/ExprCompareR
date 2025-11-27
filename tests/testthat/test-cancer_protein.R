library(testthat)

# This file includes tests for function compareCancerProtein

#-------------------------------------------------------------------------------
# Unit Tests

# Testing cancer_protein for valid inputs

test_that("compareCancerProtein warns for invalid cancer_type input", {

  # Case 1: completely invalid cancer type
  expect_error(
    compareCancerProtein("brain cancer"),
    "is not found in cancer_tissue_map\\$cancer"
  )

  # Case 2: slightly misspelled cancer type
  expect_error(
    compareCancerProtein("breastcancer"),  # missing space
    "is not found in cancer_tissue_map\\$cancer"
  )

})

#-------------------------------------------------------------------------------
# Integration Tests

# Testing if the final tibble has the necessary columns

test_that("Output table contains required columns", {

  output <- compareCancerProtein("breast cancer")
  tbl <- output$table

  expect_true(all(c(
    "ensembl", "gene", "normal_rank", "cancer_rank",
    "delta_rank", "direction"
  ) %in% names(tbl)))

})


# Testing if the function returns both a plot and vector
test_that("compareCancerProtein returns expected list structure", {

  output <- compareCancerProtein("breast cancer")

  expect_type(output, "list")
  expect_true("table" %in% names(output))
  expect_true("plot" %in% names(output))
  expect_s3_class(output$plot, "ggplot")

})

# Check that the function only keeps genes with significant expression

test_that("Filtering keeps only normal_rank between 1 and 2", {

  out <- compareCancerProtein("breast cancer")$table

  expect_true(all(out$normal_rank >= 1 & out$normal_rank <= 2))

  # This should have filtered out G1 (normal=3)
  expect_false("G1" %in% out$gene)

})

# Test if the function returns a table and plot;
# Also test if the table contains the neccessary genes

test_that("compareCancerProtein returns correct table and plot for breast cancer", {

  output <- compareCancerProtein(cancer_type = "breast cancer")

  # Check structure
  expect_true(is.list(output))
  expect_true(all(c("table", "plot") %in% names(output)))

  # Expected first 10 gene names
  expected_genes <- c(
    "TSPAN6", "DPM1", "SCYL3", "C1orf112", "NFYA",
    "NIPAL3", "ANKIB1", "BAD", "LAP3", "HECW1"
  )

  # Check table structure
  expect_s3_class(output$table, "tbl_df")
  expect_true("gene" %in% names(output$table))

  # Check that the table contains at least those 10 genes
  actual_genes <- output$table$gene[1:10]
  expect_equal(actual_genes, expected_genes)

  # Check that plot is a valid plot object
  expect_true(is.object(output$plot))
})

#-------------------------------------------------------------------------------

# [END]

