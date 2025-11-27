library(testthat)

# This file includes tests for all helper functions - get_gtex_gencode_ids, protein_expr_values, convert_to_gtex

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------

# Tests for helper function convert_to_gtex

test_that("Valid simple tissues are converted correctly", {
  input <- c("liver", "skin 1", "breast")
  expected <- tissue_map$RNA_tissue[match(tolower(input), tolower(tissue_map$protein_tissue))]

  result <- convert_to_gtex(input)
  expect_equal(result, expected)
})

test_that("Valid complex tissues are converted correctly", {
  input <- c("kidney","lung","spleen","liver","ovary","testis",
             "breast","adrenal gland","pancreas","cerebellum")
  expected <- tissue_map$RNA_tissue[match(tolower(input), tolower(tissue_map$protein_tissue))]

  result <- convert_to_gtex(input)
  expect_equal(result, expected)
})

test_that("Unknown tissues trigger warning and return NA", {
  input <- c("liver", "unknown_tissue", "skin 1")

  expect_warning(result <- convert_to_gtex(input),
                 "Some tissues could not be converted")

  # The unknown tissue should produce NA
  expect_true("unknown_tissue" %in% input[is.na(result) | result == "Unknown"])
})

test_that("Non-character input triggers error", {
  expect_error(convert_to_gtex(123), "`tissues` must be a character vector")
  expect_error(convert_to_gtex(list("liver")), "`tissues` must be a character vector")
})

test_that("Empty input returns character(0)", {
  expect_equal(convert_to_gtex(character(0)), character(0))
})

#-------------------------------------------------------------------------------
# [END]


