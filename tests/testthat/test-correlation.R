library(testthat)

# This file includes tests for function compute_correlation and correlation_genes_tissues

#-------------------------------------------------------------------------------
# Test compute_correlation

# Unit Tests

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

#-------------------------------------------------------------------------------

# Test correlation_genes_tissues

# Unit Tests

# Dummy valid inputs
valid_genes   <- c("A", "B", "C", "D", "E")
valid_tissues <- c("lung", "heart", "liver", "brain", "kidney")

test_that("Error if both gene_NAMES and tissue_NAMES are NULL", {
  expect_error(
    correlation_genes_tissues(NULL, NULL, plot_choice = "per_gene"),
    "Please provide at least five gene names AND at least five tissue names"
  )
})

test_that("Error if only gene_NAMES is provided", {
  expect_error(
    correlation_genes_tissues(valid_genes, NULL, plot_choice = "per_gene"),
    "requires both gene_NAMES and tissue_NAMES"
  )
})

test_that("Error if only tissue_NAMES is provided", {
  expect_error(
    correlation_genes_tissues(NULL, valid_tissues, plot_choice = "per_gene"),
    "requires both gene_NAMES and tissue_NAMES"
  )
})

test_that("Error if gene_NAMES has fewer than five elements", {
  expect_error(
    correlation_genes_tissues(c("A","B","C"), valid_tissues, "per_gene"),
    "gene_NAMES must contain at least five elements"
  )
})

test_that("Error if tissue_NAMES has fewer than five elements", {
  expect_error(
    correlation_genes_tissues(valid_genes, c("lung", "heart"), "per_gene"),
    "tissue_NAMES must contain at least five elements"
  )
})

test_that("Error if duplicate genes provided", {
  expect_error(
    correlation_genes_tissues(c("A","B","C","A","D"), valid_tissues, "per_gene"),
    "Duplicate gene names detected"
  )
})

test_that("Error if duplicate tissues provided", {
  expect_error(
    correlation_genes_tissues(valid_genes, c("lung","heart","lung","kidney","liver"), "per_gene"),
    "Duplicate tissue names detected"
  )
})

test_that("Error if invalid plot_choice is provided", {
  expect_error(
    correlation_genes_tissues(valid_genes, valid_tissues, plot_choice = "invalid"),
    "plot_choice is invalid"
  )
})

#-------------------------------------------------------------------------------
# Integration Tests:

# Test if compute_correlation returns a plot object for a tissue list

test_that("compute_correlation returns a plot object for per_tissue_plot", {

  tissue_list <- c("lung", "spleen", "liver", "ovary", "testis")

  output <- compute_correlation(tissue_NAMES = tissue_list)

  # Structure checks
  expect_true(is.list(output))
  expect_true("per_tissue_plot" %in% names(output))

  # Check that per_tissue_plot is a plot object
  expect_true(is.object(output$per_tissue_plot))
})

#-------------------------------------------------------------------------------

# Integration Tests:

# Test if correlation_genes_tissues returns a plot object, if plot_choice is per_gene

test_that("correlation_genes_tissues returns a plot object for per_gene plot choice", {

  tissues <- c("lung", "spleen", "liver", "ovary", "testis")
  genes <- c("MYC", "TP53", "BRCA1", "CRP", "EGFR")

  result_gene <- correlation_genes_tissues(
    gene_NAMES = genes,
    tissue_NAMES = tissues,
    plot_choice = "per_gene"
  )

  # Structure checks
  expect_true(is.list(result_gene))
  expect_true("per_gene_plot" %in% names(result_gene))

  # Check that per_gene_plot is a plot object
  expect_true(is.object(result_gene$per_gene_plot))
})

# Test if correlation_genes_tissues returns a plot object, if plot_choice is per_tissue

test_that("correlation_genes_tissues returns a plot object for per_tissue plot choice", {

  tissues <- c("lung", "spleen", "liver", "ovary", "testis")
  genes <- c("MYC", "TP53", "BRCA1", "CRP", "EGFR")

  result_tissue <- correlation_genes_tissues(
    gene_NAMES = genes,
    tissue_NAMES = tissues,
    plot_choice = "per_tissue"
  )

  # Structure checks
  expect_true(is.list(result_tissue))
  expect_true("per_tissue_plot" %in% names(result_tissue))

  # Check it's a plot object
  expect_true(is.object(result_tissue$per_tissue_plot))
})

#-------------------------------------------------------------------------------


# [END]

