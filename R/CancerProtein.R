#' Comparison of protein expression between normal and cancer tissues
#'
#' This function compares the protein expression levels between normal and
#' pathology (cancer) tissues for a specified cancer type. It calculates the
#' change in expression rank and direction ("Up", "Down", or
#' "No change") for each gene, and visualizes the shift in expression.
#'
#'@param cancer_type A character string specifying the cancer type of interest.
#' The value must be one of the valid cancer names available in
#' \code{cancer_tissue_map$cancer}. If an invalid cancer type is provided,
#' the function will stop execution and throw an error.
#'
#' @return
#' A list containing two components:
#' \itemize{
#'   \item \code{table} - A tibble of metadata with the following columns:
#'     \describe{
#'       \item{ensembl}{Ensembl gene identifier.}
#'       \item{gene}{Gene symbol.}
#'       \item{normal_rank}{Average protein expression rank in normal tissue.}
#'       \item{cancer_rank}{Weighted protein expression rank in cancer tissue.}
#'       \item{delta_rank}{Difference (\eqn{cancer\_rank - normal\_rank}).}
#'       \item{direction}{Direction of change ("Up", "Down", or "No change").}
#'     }
#'   \item \code{plot} - A \pkg{ggplot2} bar plot showing the direction and
#'   magnitude of protein rank changes for each gene.
#'   The down regulated genes are labeled in red, the up regulated genes are labeled in blue,
#'   and the genes with no change in protein expression are labeled in green.
#'
#' }
#'
#'@details
#' The function uses the \pkg{HPAanalyze} dataset (`hpaDownload`) for
#' expression data, and the internal dataset \code{\link{cancer_tissue_map}}
#' to match each cancer type to its corresponding tissue.
#'
#'@note ExprCompareR retrieves data from HPA (Thul et al., 2018) and GTEx (Lonsdale et al., 2013)
#'via packages, HPAanalyze(Tran et al., 2019) and gtexr(Warwick et al., 2025).
#'Because these functions make external API calls, some operations may take time to complete
#'typically around 30 seconds to 1 minute.
#'
#'If you are interested in obtaining the top 10 up-regulated genes, top 10 down-regulated genes or top 10 no change genes from the cancer metadata tibble, then you can filter the output$table to obtain these lists.
#'See "inst/shiny-scripts/modules/module-cancer_protein.R" lines 139-170 for help.
#'
#'
#'@examples
#'\dontrun{
#' # Example 1: Compare protein expression in a common tissue
#'
#' output <- compareCancerProtein(cancer_type = "breast cancer")
#' output$table
#' output$plot
#'}
#'
#'@references

#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
#'Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#'Springer-Verlag New York. ISBN 978-3 319-24277-4, \url{https://ggplot2.tidyverse.org.}
#'
#'Wickham H, Francois R, Henry L, Muller K, Vaughan D (2025).
#'dplyr: A Grammar of Data Manipulation. R package version 1.1.4,
#'\url{https://dplyr.tidyverse.org.}
#'
#'@import HPAanalyze
#'@import ggplot2
#'@import dplyr
#'@importFrom stats reorder
#'@importFrom stats cor median
#'@importFrom utils head
#'
#'@export
compareCancerProtein<- function(cancer_type){

  # Testing if the cancer_type exists in cancer_tissue_map$cancer
  if (!cancer_type %in% cancer_tissue_map$cancer) {
    stop("The cancer type '", cancer_type, "' is not found in cancer_tissue_map$cancer.")
  }
  # Subset pathology for protein expression related to the given cancer type
  cancer_data <- pathology %>%
    dplyr::filter(cancer == cancer_type)

  # Create a vector of gene names involved in cancer_data, so that we can query
  # the normal tissue expression of these genes
  genes_involved <- unique(cancer_data$gene)

  # Use data cancer_tissue_map to find the equivalent tissues primarily
  # involved affected by the cancer
  type_tissue <- cancer_tissue_map %>%
    dplyr::filter(cancer == cancer_type) %>%
    dplyr::pull(tissue)


  # If their is no equivalent tissue, stop analysis.
  if(length(type_tissue) == 0){
    stop("No tissue match found for ", cancer_type, ". This is an invalid cancer type.")
  }

  # Query normal tissue for genes in genes_involved and those which have tissue type, to
  # obtain normal gene expression
  normal_filtered <- normal_tissue %>%
    dplyr::filter(gene %in% genes_involved, tissue == type_tissue)

  # Converting the protein expression string values to numeric values
  normal_filtered <- normal_filtered %>%
    dplyr::mutate(rank_value = protein_expr_values(level))

  # Computing the average of different samples of a gene in the tissue
  # as we have multiple samples of a given gene.
  normal_summary <- normal_filtered %>%
    dplyr::group_by(ensembl, gene) %>%
    dplyr::summarise(normal_rank = mean(rank_value, na.rm = TRUE), .groups = "drop")

  # Compute weighted rank for cancer_data
  cancer_summary <- cancer_data %>%
    dplyr::mutate(cancer_rank = (1*high + 2*medium + 3*low + 0*`not_detected`)/
                    (high + medium + low + not_detected)) %>%
    dplyr::select(ensembl, gene, cancer_rank)


  # Compute the difference in cancer vs normal expression of all genes.
  # Assign a rank of "Up", "Down" or "No change" depending on the delta.
  # Store the metadata in a table named results.
  results <- dplyr::full_join(normal_summary, cancer_summary, by = c("ensembl", "gene")) %>%
    dplyr::mutate(
      delta_rank = cancer_rank - normal_rank,
      direction = dplyr::case_when(
        delta_rank > 0 ~ "Up",
        delta_rank < 0 ~ "Down",
        TRUE ~ "No change"
      )
    )

  # Filter genes with normal_rank between 1 and 2, to plot only genes with
  # a significant expression
  results <- results %>%
    dplyr::filter(normal_rank >= 1 & normal_rank <= 2)

  # Plot the shift in rank for all genes specifically highlight the direction.
  delta_plot <- ggplot(results, aes(x = reorder(gene, delta_rank), y = delta_rank, fill = direction)) +
    geom_col() +
    coord_flip() +
    labs(title = paste("Protein Rank Shift in", cancer_type),
         x = NULL, y = "Delta Rank (Cancer - Normal)") +
    theme_minimal() +
    theme(axis.text.y = element_blank())

  # return an object which has the results of direction and detla, and the plot of
  # delta and direction
  output <- list(table = results, plot = delta_plot)

  return (output)

}

# [END]
