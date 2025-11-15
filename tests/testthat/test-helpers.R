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


