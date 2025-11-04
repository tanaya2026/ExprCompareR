library(testthat)

test_that("compute_correlation handles invalid inputs with warnings", {

  # Case 1: both NULL
  expect_warning(
    compute_correlation(gene_NAMES = NULL, tissue_NAMES = NULL),
    "Please provide at least five gene names or tissue names"
  )

  # Case 2: both provided but too few elements
  expect_warning(
    compute_correlation(gene_NAMES = c("BRCA1", "TP53"),
                        tissue_NAMES = c("liver", "lung")),
    "Please provide at least five gene names or tissue names"
  )

  # Case 3: only gene_NAMES provided but <5 elements
  expect_warning(
    compute_correlation(gene_NAMES = c("BRCA1", "TP53", "MYC")),
    "Need at least five gene names"
  )

  # Case 4: only tissue_NAMES provided but <5 elements
  expect_warning(
    compute_correlation(tissue_NAMES = c("liver", "lung", "spleen")),
    "Need at least five tissue names"
  )

})
