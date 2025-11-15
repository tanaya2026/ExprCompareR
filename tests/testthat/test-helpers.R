library(testthat)

# Tests for helper function get_gtex_gencode_ids

test_that("get_gtex_gencode_ids works for valid input", {
  # Pick 1–2 known valid genes from your gene_symbols_list
  valid_genes <- head(gene_symbols_list, 2)

  # Mock GTEx response if possible, otherwise this test calls GTEx
  gencode_ids <- get_gtex_gencode_ids(valid_genes)

  expect_type(gencode_ids, "character")
  expect_length(gencode_ids, length(valid_genes))
})

test_that("Non-character input throws error", {
  expect_error(get_gtex_gencode_ids(123),
               "`gene_symbols` must be a character vector")
  expect_error(get_gtex_gencode_ids(list("MYC")),
               "`gene_symbols` must be a character vector")
})

test_that("Invalid gene symbols throw error", {
  expect_error(get_gtex_gencode_ids(c("INVALIDGENE")),
               "All `gene_symbols` must be valid symbols listed in `gene_symbols_list`")
})

test_that("Empty input returns empty vector", {
  expect_equal(get_gtex_gencode_ids(character(0)), character(0))
})

#________________________________________________________________________________

# Tests for helper function protein_expr_values

test_that("Valid protein expression levels are converted correctly", {
  input <- c("Not detected", "n/a", "na", "NA", "Low", "Medium", "High", NA)
  expected <- c(0, 0, 0, 0, 1, 2, 3, NA)

  result <- protein_expr_values(input)
  expect_equal(result, expected)
})

test_that("Non-character input triggers an error", {
  expect_error(protein_expr_values(123), "`protein_expression` must be a character vector")
  expect_error(protein_expr_values(list("Low")), "`protein_expression` must be a character vector")
})

test_that("Unknown string values trigger an error", {
  expect_error(protein_expr_values(c("Very High", "Low")),
               "Invalid protein expression values found: Very High")
})

test_that("Empty input returns numeric(0)", {
  expect_equal(protein_expr_values(character(0)), numeric(0))
})

test_that("NA values are preserved", {
  input <- c(NA, "Low", NA)
  result <- protein_expr_values(input)
  expect_equal(result, c(NA, 1, NA))
})

test_that("Function handles vector of only recognized zeros correctly", {
  input <- c("Not detected", "n/a", "na", "NA")
  expect_equal(protein_expr_values(input), c(0, 0, 0, 0))
})

test_that("Function handles vector of only High, Medium, Low correctly", {
  input <- c("Low", "Medium", "High")
  expect_equal(protein_expr_values(input), c(1, 2, 3))
})


