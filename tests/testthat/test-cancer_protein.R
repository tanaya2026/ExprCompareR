library(testthat)


#_______________________________________________________________________________

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


#_______________________________________________________________________________

# Testing functionality of cancer_protein

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

#_______________________________________________________________________________


