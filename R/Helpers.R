#' Converts Gene Symbols to GENCODE IDs
#'
#' The function given a vector of gene symbols, converts them into their
#' respective GENCODE IDs. This helper function is useful when querying RNA expression
#' data of a given gene, using the \pkg{gtexr} package, which requires input genes
#' to be provided as GENCODE IDs.
#' It serves as an internal utility within the ExprCompareR package, allowing
#' users to easily map gene symbols (from the Human Protein Atlas) to their
#' GENCODE equivalents required for GTEx queries.
#'
#'@param gene_symbols A character vector of gene symbols to query.
#' Valid symbols are listed in \code{\link{gene_symbols_list}}.
#' These genes typically have protein expression data but may not always
#' have corresponding RNA expression data. Missing or unmatched genes will
#' trigger a warning and be excluded from the output.
#'
#'@return A character vector of corresponding GENCODE IDs.
#'Genes not found in the GTEx database are omitted with a warning.
#'
#'@details
#'The function retrieves metadata from the GTEx database via the
#'\code{gtexr::get_genes()} function and extracts the \code{gencodeId} values
#'corresponding to the provided gene symbols.
#'
#'@examples
#'\dontrun{
#' # Example 1: Retrieve GENCODE IDs for common genes.
#'
#' my_genes <- c("MYC", "TP53", "BRCA1")
#'
#' gencode_ids <- get_gtex_gencode_ids(gene_symbols = my_genes)
#' gencode_ids
#'
#'
#' # Example 2: Include a gene symbol which does not exist (will issue a warning).
#'
#' my_genes <- c("CRP", "EGFR", "UNKNOWN")
#'
#' gencode_ids <- get_gtex_gencode_ids(gene_symbols = my_genes)
#' gencode_ids
#'}
#'
#' @references
#'
#'Lonsdale, J., Thomas, J., Salvatore, M. et al.
#'The Genotype-Tissue Expression (GTEx) project. Nat Genet 45, 580-585 (2013).
#'https://doi.org/10.1038/ng.2653
#'
#'Thul PJ, Lindskog C. The human protein atlas: A spatial map of the human proteome.
#'Protein Sci. 2018 Jan;27(1):233-244. doi: 10.1002/pro.3307.
#'Epub 2017 Oct 10. PMID: 28940711; PMCID: PMC5734309.
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'https://doi.org/10.1186/s12859-019-3059-z
#'
#'Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). "gtexr: A
#'convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API."
#'Journal of Open Source So ware, 10(109), 8249. ISSN
#'2475-9066, doi:10.21105/joss.08249, gigs v0.2.1.
#'
#'@import HPAanalyze
#'@import gtexr
#'@import dplyr
#'@importFrom stats cor median
#'@importFrom utils head


get_gtex_gencode_ids <- function(gene_symbols) {

  # Testing Input validation
  # Test if all inputs are characters
  if (!is.character(gene_symbols)) {
    stop("`gene_symbols` must be a character vector of gene symbols.")
  }

  # Test if all gene_symbols are valid i.e. in gene_symbols_list
  if (!all(gene_symbols %in% gene_symbols_list)) {
    stop("All `gene_symbols` must be valid symbols listed in `gene_symbols_list`.")
  }

  # If the user inputs an empty vector, then return an empty character vector
  if (length(gene_symbols) == 0) {
    return(character(0))
  }



  # Call get_genes() from gtexr with user-provided gene symbols,
  # to get gene metadata
  gene_info <- suppressMessages(gtexr::get_genes(gene_symbols))

  # Extract the IDs from the gene metadata, and create a vector of Gencode IDs
  gencode_vector <- sapply(gene_symbols, function(gene) {
    match_row <- gene_info %>% filter(geneSymbol == gene)
    # If the given gene's GENCODE ID does not exist, return NA
    if (nrow(match_row) == 0) {
      warning("Gene not found in GTEx: ", gene)
      return(NA_character_)
    } else {
      return(match_row$gencodeId[1])
    }
  })

  # Remove symbols of genes from vector
  gencode_vector <- unname(gencode_vector)

  # Remove NA values, as these are invalid genes
  gencode_vector <- gencode_vector[!is.na(gencode_vector)]

  # Return the vector of GENCODE IDs
  return(gencode_vector)
}

#' Converts HPA protein expression levels to numeric values
#'
#' The Human Protein Atlas (HPA) represents protein expression levels as strings
#' such as "Not detected", "Low", "Medium", or "High". This helper function converts
#' these string labels into numeric values so that they can be used in downstream
#' quantitative analyses, such as comparing RNA and protein expression.
#'
#'@param protein_expression A character vector of protein expression levels,
#'   as returned by HPAanalyze. Values can include "Not detected", "Low", "Medium",
#'   "High", or NA.
#'
#'@return A numeric vector of the same length as \code{protein_expression}
#'with the following mapping:
#' \itemize{
#'   \item "Not detected", "n/a", "na", "NA" -> 0
#'   \item "Low" -> 1
#'   \item "Medium" -> 2
#'   \item "High" -> 3
#'   \item NA or unrecognized values -> NA
#' }
#'
#'@examples
#'\dontrun{
#' # Example 1: Retrieve numeric values for HPA expression values
#'
#' protein_levels <- c("High", "Medium", "Low", "Not detected", NA)
#'
#' numeric_values <- protein_expr_values(protein_expression = protein_levels)
#' numeric_values
#'}
#'
#' @references
#'
#'Thul PJ, Lindskog C. The human protein atlas: A spatial map of the human proteome.
#'Protein Sci. 2018 Jan;27(1):233-244. doi: 10.1002/pro.3307.
#'Epub 2017 Oct 10. PMID: 28940711; PMCID: PMC5734309.
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'https://doi.org/10.1186/s12859-019-3059-z

protein_expr_values <- function(protein_expression){

  # Generate the final result vector
  result <- numeric(length(protein_expression))

  # Loop over each element, and convert the string to a numeric value
  for (i in seq_along(protein_expression)){
    v<-protein_expression[i]

    if (is.na(v)) {
      result[i] <- NA
    }
    else if (v %in% c("Not detected", "n/a", "na", "NA")) {
      result[i]<- 0
    }
    else if(v=="Low") {
      result[i]<- 1
    }
    else if(v=="Medium") {
      result[i]<- 2
    }
    else if(v=="High") {
      result[i]<- 3
    }
    else {
      result[i] = NA
    }

  }
  # Return the final vector with numeric values
  return (result)
}


#'Convert HPA tissue names to GTEx tissue names
#'
#'Converts a vector of tissue names from the Human Protein Atlas (HPA) format
#'to GTEx-compatible tissue names using the \code{tissue_map} dataset. The tissue
#'names that each package and database use are different, and this helper function
#'allows easy conversion for analysis.
#'If a tissue cannot be converted, it returns "Unknown" and issues a warning.
#'
#'@param tissues Character vector of HPA tissue names.
#'Valid tissue names are listed in \code{tissue_map$protein_tissue}.
#'Any name not in this list will return  and issue a warning.
#'
#'@return A character vector of GTEx tissue names corresponding to input tissues.
#
#'
#'@examples
#'\dontrun{
#' # Example 1: Convert list of simple tissues
#'
#' example_input <- c("liver", "skin 1", "breast")
#' gtexr_tissues <- convert_to_gtex(tissues = example_input)
#' gtexr_tissues
#'
#' # Example 2: Convert list of complex tissues
#'
#' example_input <- c(
#' "kidney","lung", "spleen", "liver", "ovary", "testis",
#' "breast","adrenal gland", "pancreas", "cerebellum")
#' gtexr_tissues <- convert_to_gtex(tissues = example_input)
#' gtexr_tissues
#'
#'}
#'
#'@references
#'
#'Lonsdale, J., Thomas, J., Salvatore, M. et al.
#'The Genotype-Tissue Expression (GTEx) project. Nat Genet 45, 580-585 (2013).
#'https://doi.org/10.1038/ng.2653
#'
#'Thul PJ, Lindskog C. The human protein atlas: A spatial map of the human proteome.
#'Protein Sci. 2018 Jan;27(1):233-244. doi: 10.1002/pro.3307.
#'Epub 2017 Oct 10. PMID: 28940711; PMCID: PMC5734309.
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'https://doi.org/10.1186/s12859-019-3059-z
#'
#'Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). "gtexr: A
#'convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API."
#'Journal of Open Source So ware, 10(109), 8249. ISSN
#'2475-9066, doi:10.21105/joss.08249, gigs v0.2.1.


convert_to_gtex <- function(tissues) {
    # Lookup in tissue_map dataframe
    idx <- match(tolower(tissues), tolower(tissue_map$protein_tissue))
    converted <- tissue_map$RNA_tissue[idx]

    # Warn if any tissue failed
    if (any(is.na(converted))) {
      warning("Some tissues could not be converted: ",
              paste(tissues[is.na(converted)], collapse = ", "))
    }

    return(unname(converted))
  }




# [END]
