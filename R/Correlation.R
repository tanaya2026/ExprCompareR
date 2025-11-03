#' Converts protein expression data queried from HPAanalyze, to numeric values
#'
#' A function that converts the protein expression data of a given gene/tissue
#' from HPAanalyze, into numeric values. HPAnalyze returns protein expression as
#' strings, for e.g.("low", "medium", etc). In order to determine the pearson
#' correlation between the RNA expression and protein expression, we convert these
#' strings into a numeric value.
#'



#________________________________________________________________________________


compute_correlation <- function(gene_NAMES = NULL, tissue_NAMES = NULL){

  # Case 1: Inputs, gene_NAMES and tissue_NAMES are both NULL
  if (is.null(gene_NAMES) && is.null(tissue_NAMES)){
    warning("Please provide at least five gene names or tissue names")
    return (NA)
  }


  # Case 2: Not enough entries in gene_NAMES and tissue_NAMES
  if ((!is.null(gene_NAMES) && length(gene_NAMES) < 5) &&
      (!is.null(tissue_NAMES) && length(tissue_NAMES) < 5)) {
    warning("Please provide at least five gene names or tissue names.")
    return(NA)
  }

  # Case 3: if only one list is provided, ensure it has >= 3 elements
  if (!is.null(gene_NAMES) && is.null(tissue_NAMES) && length(gene_NAMES) < 5) {
    warning("Need at least five gene names.")
    return(NA)
  }

  if (is.null(gene_NAMES) && !is.null(tissue_NAMES) && length(tissue_NAMES) < 5) {
    warning("Need at least five tissue names.")
    return(NA)
  }

  # Main_logic:

  # Needs to give gene IDs -> bioMART

  if (is.null(tissue_NAMES) || length(tissue_NAMES) == 0){
    RNA_expr_results <- list()
  }
}


#_________________________________________________________________________

correlation_genes_only <- function(gene_list){






}
