library(testthat)


test_that("compute_correlation handles all 7 invalid input cases", {

  # Case 1: both gene_NAMES and tissue_NAMES are NULL
  expect_error(
    compute_correlation(gene_NAMES = NULL, tissue_NAMES = NULL),
    "Please provide at least five gene names OR tissue names"
  )

  # Case 2: both provided but too few elements in both
  expect_error(
    compute_correlation(
      gene_NAMES = c("BRCA1", "TP53"),
      tissue_NAMES = c("liver", "lung")
    ),
    "Please provide at least five gene names OR tissue names."
  )

  # Case 3: only gene_NAMES provided but <5 elements
  expect_error(
    compute_correlation(gene_NAMES = c("BRCA1", "TP53", "MYC")),
    "Need at least five gene names."
  )

  # Case 4: only tissue_NAMES provided but <5 elements
  expect_error(
    compute_correlation(tissue_NAMES = c("liver", "lung", "spleen")),
    "Need at least five tissue names."
  )

  # Case 5: both lists provided (function supports only one list)
  expect_error(
    compute_correlation(
      gene_NAMES = c("BRCA1", "TP53", "MYC", "CRP", "EGFR"),
      tissue_NAMES = c("liver", "lung", "spleen", "ovary", "testis")
    ),
    "This function only supports either genes OR tissues, not both"
  )

  # Case 6: duplicate gene names
  expect_error(
    compute_correlation(gene_NAMES = c("BRCA1", "TP53", "BRCA1", "MYC", "EGFR")),
    "Duplicate gene names detected: BRCA1"
  )

  # Case 7: duplicate tissue names
  expect_error(
    compute_correlation(tissue_NAMES = c("liver", "lung", "liver", "spleen", "ovary")),
    "Duplicate tissue names detected: liver"
  )

})

