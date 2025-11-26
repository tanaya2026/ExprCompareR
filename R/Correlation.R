#' Compute Correlation for Genes or Tissues
#'
#'
#' This function serves as a wrapper to compute and visualize the Spearman correlation
#' between RNA and protein expression between a user's gene list **OR** tissue list of interest and plots the results.
#' It requires at least five entries in the provided list (`gene_NAMES` or `tissue_NAMES`)
#' to perform the computation. Depending on the input, it calls either
#' `correlation_genes_only()` or `correlation_tissues_only()` internally.
#'
#' @note
#' This function **ONLY** computes correlations either for a list of genes **OR** a list of tissues.
#' If you wish to input a list of genes **AND** list of tissues, utilize function `correlation_genes_tissues`, or
#' see `help` documentation `?correlation_genes_tissues`
#'
#' @param gene_NAMES A character vector of gene names.
#' Optional; defaults to \code{NULL}.
#' If provided, must contain at least five elements.
#' Valid gene symbols for the \code{gene_NAMES} argument can be accessed from
#' \code{gene_symbols_list}.
#'
#' @param tissue_NAMES A character vector of tissue names.
#' Optional; defaults to \code{NULL}.
#' If provided, must contain at least five elements.
#' Valid tissue names as for the \code{tissue_NAMES} argument can be accessed from
#' \code{tissue_map$protein_tissue}
#'
#' @return A plot object:
#'   - If `gene_NAMES` is provided and `tissue_NAMES` is NULL, returns the `per_gene_plot`  of Spearman correlations from `correlation_genes_only()`. The Spearman correlation values for each gene is labeled in the plot.
#'   - If `tissue_NAMES` is provided and `gene_NAMES` is NULL, returns the `per_tissue_plot` of Spearman correlations from `correlation_tissues_only()`. The Spearman correlation values for each tissue is labeled in the plot.
#'   - Stops execution and throws an error if input criteria is not met.
#'
#' @details
#' - This function utilizes data/tissue_list_RNA.rda and data/tissue_list_protein.rda
#'
#' The function validates the input before computing correlations:
#' - If both `gene_NAMES` and `tissue_NAMES` are NULL, it stops execution and throws an error
#' - If one list is provided but has fewer than 5 elements, it stops execution and throws an error
#' - If both lists are provided with 5>= elements each, it stops execution and throws an error
#' - If duplicates exist in either list, it stops execution and throws an error
#' - If only one list is provided, it must contain at least 5 elements.
#' - Plot key - Spearman correlation values are represented by color intensity. High positive correlations are indicated in dark blue,
#' low correlations around zero are white, and high negative correlations are indicated in dark red.
#' The shade intensity reflects the magnitude of the correlation.
#'
#'
#' @note ExprCompareR retrieves data from HPA (Thul et al., 2018) and GTEx (Lonsdale et al., 2013)
#'via packages, HPAanalyze(Tran et al., 2019) and gtexr(Warwick et al., 2025).
#'Because these functions make external API calls, some operations may take time to complete
#'typically around 30 seconds to 1 minute.
#'
#'
#' @examples
#' \dontrun{
#' # Example 1: Using gene names only
#'
#' gene_list <- c("MYC", "TP53", "BRCA1", "CRP", "EGFR")
#' output<- compute_correlation(gene_NAMES = gene_list)
#' output$per_gene_plot
#'
#' # Example 2: Using tissue names only
#'
#' tissue_list = c("lung", "spleen", "liver", "ovary", "testis")
#' output<- compute_correlation(tissue_NAMES = tissue_list)
#' output$per_tissue_plot
#'
#'
#' # Example 3: Invalid input (too few elements)
#' compute_correlation(gene_NAMES = c("BRCA1", "TP53"))
#' }
#'
#'@references
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
#'Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). "gtexr: A
#'convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API."
#'Journal of Open Source Software, 10(109), 8249. ISSN 2475-9066,
#'\url{doi:10.21105/joss.08249}, gigs v0.2.1.
#'
#'Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#'Springer-Verlag New York. ISBN 978-3 319-24277-4, \url{https://ggplot2.tidyverse.org.}
#'
#'Wickham H, Francois R, Henry L, Muller K, Vaughan D (2025).
#'dplyr: A Grammar of Data Manipulation. R package version 1.1.4,
#'\url{https://dplyr.tidyverse.org.}
#'
#'Wickham H, Henry L (2025). purrr: Functional Programming Tools.
#'R package version 1.1.0, \url{https://purrr.tidyverse.org/.}
#'
#'Wickham H, Vaughan D, Girlich M (2025). tidyr: Tidy Messy Data.
#'R package version 1.3.1, \url{https://tidyr.tidyverse.org.}
#'
#'
#'@importFrom stats cor median
#'@importFrom utils head
#' @export

compute_correlation <- function(gene_NAMES = NULL, tissue_NAMES = NULL){

  # Case 1: Inputs, gene_NAMES and tissue_NAMES are both NULL
  if (is.null(gene_NAMES) && is.null(tissue_NAMES)){
    stop("Please provide at least five gene names OR tissue names")
  }

  # Case 2: Not enough entries in gene_NAMES or tissue_NAMES
  if ((!is.null(gene_NAMES) && length(gene_NAMES) < 5) &&
      (!is.null(tissue_NAMES) && length(tissue_NAMES) < 5)) {
    stop("Please provide at least five gene names OR tissue names.")
  }

  # Case 3: If only one list is provided, ensure it has >= 5 elements
  if (!is.null(gene_NAMES) && is.null(tissue_NAMES) && length(gene_NAMES) < 5) {
    stop("Need at least five gene names.")
  }

  # Case 4: If only one list is provided, ensure it has >= 5 elements
  if (is.null(gene_NAMES) && !is.null(tissue_NAMES) && length(tissue_NAMES) < 5) {
    stop("Need at least five tissue names.")
  }

  # Case 5: If both lists are provided, then stop
  if (!is.null(gene_NAMES) && !is.null(tissue_NAMES)) {
  stop("This function only supports either genes OR tissues, not both. Use `correlation_genes_tissues` for both.")
  }

  # Case 6: Check for duplicates in gene_NAMES
  if (!is.null(gene_NAMES) && any(duplicated(gene_NAMES))) {
    dupes <- unique(gene_NAMES[duplicated(gene_NAMES)])
    stop(
      paste0(
        "Duplicate gene names detected: ",
        paste(dupes, collapse = ", "),
        ".\nPlease remove duplicates before proceeding."
      )
    )
  }

  # Case 7: Check for duplicates in tissue_NAMES
  if (!is.null(tissue_NAMES) && any(duplicated(tissue_NAMES))) {
    dupes <- unique(tissue_NAMES[duplicated(tissue_NAMES)])
    stop(
      paste0(
        "Duplicate tissue names detected: ",
        paste(dupes, collapse = ", "),
        ".\nPlease remove duplicates before proceeding."
      )
    )
  }



  # If only gene_NAMES are provided, call function correlation_genes_only
  if (is.null(tissue_NAMES) && length(gene_NAMES)>=5){
    output <-correlation_genes_only(gene_list = gene_NAMES)
    return(output)
  }

  # If only tissue_NAMES are provided, call function correlation_tissues_only
  else if (is.null(gene_NAMES) && length(tissue_NAMES)>=5){
    output <- correlation_tissues_only(tissue_NAMES)
    return(output)
  }
}




#' Compute Spearman Correlation Between RNA and Protein Expression for Genes
#'
#' This function calculates the Spearman correlation between RNA and protein expression levels
#' for a given list of genes across multiple standard tissues. It also generates a bar plot
#' visualizing the correlation values for each gene.
#'
#' @param gene_list A character vector of gene symbols. Only valid genes will be processed.
#' Valid gene symbols for the \code{gene_NAMES} argument can be accessed from
#' \code{gene_symbols_list}.
#'
#' @return A list containing:
#' \describe{
#'   \item{per_gene_plot}{A \code{ggplot2} object showing Spearman correlation values
#'   for each gene across tissues.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Converts gene symbols to GENCODE IDs using \code{get_gtex_gencode_ids}.
#'   \item Queries RNA expression data for standard tissues \code{tissue_list_RNA} using \code{gtexr::get_gene_expression}.
#'   \item Aggregates RNA expression by averaging multiple samples per gene per tissue.
#'   \item Downloads and filters normal tissue protein expression data using \code{HPAanalyze::hpaDownload}.
#'   \item Checks for missing gene-tissue combinations in protein data, issuing a warning and filling missing entries with NA.
#'   \item Converts protein expression levels to numeric values using \code{protein_expr_values}.
#'   \item Aligns RNA and protein data and creates a separate tibble for each gene.
#'   \item Computes Spearman correlations between RNA and protein expression across tissues.
#'   \item Generates a bar plot displaying the Spearman correlation for each gene.
#'
#' - Plot key - Spearman correlation values are represented by color intensity. High positive correlations are indicated in dark blue,
#' low correlations around zero are white, and high negative correlations are indicated in dark red.
#' The shade intensity reflects the magnitude of the correlation.
#' }
#'
#'@note ExprCompareR retrieves data from HPA (Thul et al., 2018) and GTEx (Lonsdale et al., 2013)
#'via packages, HPAanalyze(Tran et al., 2019) and gtexr(Warwick et al., 2025).
#'Because these functions make external API calls, some operations may take time to complete
#'typically around 30 seconds to 1 minute.
#'
#' @examples
#' \dontrun{
#'
#' # Example 1: Using valid gene names only
#'
#' gene_list <- c("MYC", "TP53", "BRCA1", "CRP", "EGFR")
#' output<- correlation_genes_only(gene_NAMES = gene_list)
#' output$per_gene_plot
#'
#' }
#'
#' @references
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
#'Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). "gtexr: A
#'convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API."
#'Journal of Open Source Software, 10(109), 8249. ISSN 2475-9066,
#'\url{doi:10.21105/joss.08249}, gigs v0.2.1.
#'
#'Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#'Springer-Verlag New York. ISBN 978-3 319-24277-4, \url{https://ggplot2.tidyverse.org.}
#'
#'Wickham H, Francois R, Henry L, Muller K, Vaughan D (2025).
#'dplyr: A Grammar of Data Manipulation. R package version 1.1.4,
#'\url{https://dplyr.tidyverse.org.}
#'
#'Wickham H, Henry L (2025). purrr: Functional Programming Tools.
#'R package version 1.1.0, \url{https://purrr.tidyverse.org.}
#'
#'Wickham H, Vaughan D, Girlich M (2025). tidyr: Tidy Messy Data.
#'R package version 1.3.1, \url{https://tidyr.tidyverse.org.}
#'
#' @import HPAanalyze
#' @import gtexr
#' @import ggplot2
#' @import tidyr
#' @import purrr
#' @import dplyr
#' @importFrom stats cor median
#' @importFrom utils head
#'
#' @keywords internal

correlation_genes_only <- function(gene_list){
  # Convert the gene_list to GENCODE IDs using helper function get_gtex_gencode_ids
  gene_ids <- get_gtex_gencode_ids(gene_list)


  # Querying gtexr to obtain RNA expression for the standard RNA tissue list
  RNA_expr<- suppressMessages(gtexr::get_gene_expression(gencodeIds = gene_ids, tissueSiteDetailIds = tissue_list_RNA))


  # Computing the average of different samples of a gene in the tissue
  # as we have multiple samples of a given gene.
  RNA_expr <- RNA_expr %>%
    mutate(avg_expr = map_dbl(data, ~ mean(.x, na.rm = TRUE)))

  # Filter protein data to only include genes from gene_list and tissues in
  # tissue_list_protein
  protein_levels <- normal_tissue %>%
    dplyr::filter(gene %in% gene_list, tissue %in% tissue_list_protein) %>%
    dplyr::select(gene, tissue, level)

  # As there are multiple samples, we summarize counts per gene-tissue and
  # suppress output
  invisible(protein_levels %>%
              dplyr::group_by(gene, tissue) %>%
              dplyr::summarise(n = n(), .groups = "drop"))

  # Check if there is a gene not found in a certain tissue

  # Create all possible gene-tissue combinations
  all_combos <- expand.grid(gene = gene_list, tissue = tissue_list_protein, stringsAsFactors = FALSE)

  # Identify missing combinations
  missing_combos <- dplyr::anti_join(all_combos, protein_levels, by = c("gene", "tissue"))

  # Check for missing combos and issue a warning if any exist
  if (nrow(missing_combos) != 0) {
    warning(
      paste(
        "Some gene-tissue combinations are missing in protein data.",
        "These missing combos indicate that specific genes of interest are not present in the",
        "standard tissue list.", "This may affect your final plot.",
        "You can either remove these genes from your analysis, replace these values or they
        will be filled by NA values.",
        sep = "\n"
      )
    )

    # Change the missing genes values to NA
    missing_combos <- missing_combos %>%
      dplyr::mutate(level = NA)

    # Append missing rows to protein_levels
    protein_levels <- dplyr::bind_rows(protein_levels, missing_combos)
  }

  # Converting the protein expression string values to numeric values
  protein_levels<- protein_levels %>%
    mutate(level = protein_expr_values(level))

  # Renaming columns for consistency and ease of use
  RNA_expr_clean <- RNA_expr %>%
    rename(tissue = tissueSiteDetailId, gene = geneSymbol, RNA = avg_expr)

  protein_levels_clean <- protein_levels %>%
    rename(protein = level)

  # Aggregate protein by gene and tissue (mean)
  protein_summary <- protein_levels_clean %>%
    dplyr::group_by(gene, tissue) %>%
    dplyr::summarise(protein = mean(protein, na.rm = TRUE), .groups = "drop")

  # Convert RNA tissues to protein_tissue
  RNA_expr_clean <- RNA_expr_clean %>%
    dplyr::left_join(tissue_map, by = c("tissue" = "RNA_tissue")) %>%
    dplyr::select(-tissue) %>%
    dplyr::rename(tissue = protein_tissue)

  # Create a tibble for each gene, in which every row is the tissue from tissue_list_RNA.
  # Each tibble has two columns, RNA and protein, and each value represents that gene's
  # RNA and protein expression in that tissue.

  gene_tables <- purrr::map(unique(RNA_expr_clean$gene), function(g) {

    RNA_sub <- RNA_expr_clean %>%
      dplyr::filter(gene == g) %>%
      dplyr::select(tissue, RNA)

    protein_sub <- protein_summary %>%
      dplyr::filter(gene == g) %>%
      dplyr::select(tissue, protein)

    final_table <- dplyr::inner_join(RNA_sub, protein_sub, by = "tissue")

    return(final_table)
  })

  names(gene_tables) <- unique(RNA_expr_clean$gene)


  # Calculating Spearman Correlation for each tibble in df, in order to compare
  # the correlation between RNA and protein expression across tissues.

  gene_correlations <- purrr::map_dbl(gene_tables, function(final_table) {
    cor(final_table$RNA, final_table$protein, use = "complete.obs", method = "spearman")
  })

  # Convert to a tibble for easy viewing
  gene_correlations_tbl <- tibble::tibble(
    gene = names(gene_correlations),
    spear_corr = gene_correlations
  )


  # Generating a plot using ggplot2 to visualize the Spearman correlation of each gene.
  final_plot <- ggplot(gene_correlations_tbl, aes(x = gene, y = spear_corr, fill = spear_corr)) +
    geom_col() +
    geom_text(aes(label = round(spear_corr, 2)), vjust = -0.5) +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0,   name = "Spearman Correlation\n(red = high negative, white = ~0, blue = high positive)\n(Color intensity reflects magnitude)"
) +
    labs(
      title = "Spearman Correlation between RNA and Protein Expression \n of genes of interest across various tissues",
      x = "Genes",
      y = "Spearman Correlation Value",
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )

  # Returning an object output which contains the final plot we created
  output <- list(per_gene_plot = final_plot)

  return (output)

}






#' Compute Spearman Correlation Between RNA and Protein Expression for Tissues
#'
#' This function calculates the Spearman correlation between RNA and protein expression levels
#' for the top expressed genes in each tissue provided. It generates a bar plot visualizing
#' the correlation for each tissue.
#'
#' @param tissue_NAMES A character vector of tissue names. Only valid tissues (matching
#'   the standard protein tissue list) will be processed.
#' Valid tissue names as for the \code{tissue_NAMES} argument can be accessed from
#' \code{tissue_map$protein_tissue}
#'
#' @return A list containing:
#' \describe{
#'   \item{per_tissue_plot}{A \code{ggplot2} object showing Spearman correlation values
#'   for each tissue.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Downloads and subsets normal tissue protein expression data using \code{HPAanalyze::hpaDownload}.
#'   \item Queries RNA expression data from GTEx for the top expressed genes in each tissue.
#'   \item Filters protein expression data to include only the top expressed genes.
#'   \item Converts protein expression levels to numeric values using \code{protein_expr_values}.
#'   \item Aggregates protein expression by median per gene per tissue.
#'   \item Creates a combined table of RNA and protein expression per tissue.
#'   \item Computes Spearman correlations between RNA and protein expression for each tissue.
#'   \item Generates a bar plot displaying Spearman correlations across tissues.
#'
#'
#' - Plot key - Spearman correlation values are represented by color intensity. High positive correlations are indicated in dark blue,
#' low correlations around zero are white, and high negative correlations are indicated in dark red.
#' The shade intensity reflects the magnitude of the correlation.
#'
#' }
#'
#'@note ExprCompareR retrieves data from HPA (Thul et al., 2018) and GTEx (Lonsdale et al., 2013)
#'via packages, HPAanalyze(Tran et al., 2019) and gtexr(Warwick et al., 2025).
#'Because these functions make external API calls, some operations may take time to complete
#'typically around 30 seconds to 1 minute.
#'
#' @examples
#' \dontrun{
#' # Example 1: Using tissue names only
#'
#' tissue_list = c("lung", "spleen", "liver", "ovary", "testis")
#' output<- correlation_tissues_only(tissue_NAMES = tissue_list)
#' output$per_tissue_plot
#' }
#'
#'@references
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
#'Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). "gtexr: A
#'convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API."
#'Journal of Open Source Software, 10(109), 8249. ISSN 2475-9066,
#'\url{doi:10.21105/joss.08249}, gigs v0.2.1.
#'
#'Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#'Springer-Verlag New York. ISBN 978-3 319-24277-4, \url{https://ggplot2.tidyverse.org.}
#'
#'Wickham H, Francois R, Henry L, Muller K, Vaughan D (2025).
#'dplyr: A Grammar of Data Manipulation. R package version 1.1.4,
#'\url{https://dplyr.tidyverse.org.}
#'
#'Wickham H, Henry L (2025). purrr: Functional Programming Tools.
#'R package version 1.1.0, \url{https://purrr.tidyverse.org.}
#'
#'Wickham H, Vaughan D, Girlich M (2025). tidyr: Tidy Messy Data.
#'R package version 1.3.1, \url{https://tidyr.tidyverse.org.}
#'
#'
#' @import HPAanalyze
#' @import gtexr
#' @import ggplot2
#' @import tidyr
#' @import purrr
#' @import dplyr
#' @importFrom stats cor median
#' @importFrom utils head

#' @keywords internal

correlation_tissues_only <- function(tissue_NAMES){

  # Creating final lists to store RNA and protein expression data
  RNA_expr_list <- list()
  Protein_expr_list <- list()

  # Querying for RNA expression data for each tissue
  for (tissue in tissue_NAMES){

    current_RNA_tissue = convert_to_gtex(tissue)

    RNA_results_tissue <- suppressWarnings( suppressMessages(gtexr::get_top_expressed_genes(
      tissueSiteDetailId = current_RNA_tissue,
      itemsPerPage = 200)
    )
    )

    # Extracting the gene names in the top expressed genes of a given tissue
    top_gene_names <- head(RNA_results_tissue$geneSymbol, 200)

    # Save RNA data as a tibble
    RNA_expr_list[[current_RNA_tissue]] <- RNA_results_tissue %>%
      dplyr::select(geneSymbol, median) %>%
      arrange(desc(median)) %>%
      as_tibble()


    protein_tissue <- tissue

    # Query protein expression data with the given gene names and tissue
    protein_tbl <- normal_tissue %>%
      filter(tissue == !!protein_tissue, gene %in% top_gene_names) %>%
      dplyr::select(gene, tissue, level) %>%
      as_tibble()

    # Save protein expression data as a tibble
    Protein_expr_list[[protein_tissue]] <- protein_tbl


  }


  # Converting the protein expression string values to numeric values
  Protein_expr_list <- lapply(Protein_expr_list, function(tbl) {
    tbl %>%
      dplyr::mutate(level = protein_expr_values(level))
  })


  # As there are multiple samples, we summarize counts per gene-tissue
  Protein_expr_median <- lapply(Protein_expr_list, function(tbl) {
    tbl %>%
      dplyr::group_by(gene, tissue) %>%
      dplyr::summarise(
        median_level = median(level, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
  })


  # Create a table for RNA and protein expression, if protein does not exist,
  # ignore the protein data.

  # Create final tibbles with RNA and protein per tissue
  Final_expr_list <- lapply(tissue_NAMES, function(tissue) {

    # Convert tissue to GTEx name for RNA
    RNA_tissue <- convert_to_gtex(tissue)

    # Get RNA tibble (geneSymbol, median)
    RNA_tbl <- RNA_expr_list[[RNA_tissue]] %>%
      dplyr::rename(RNA = median)

    # Get protein tibble (median per gene)
    protein_tbl <- Protein_expr_median[[tissue]] %>%
      dplyr::rename(protein = median_level)

    # Merge by gene
    merged_tbl <- dplyr::inner_join(RNA_tbl, protein_tbl, by = c("geneSymbol" = "gene")) %>%
      dplyr::select(geneSymbol, RNA, protein)

    return(merged_tbl)
  })

  # Name the list by tissues
  names(Final_expr_list) <- tissue_NAMES


  # Calculate the Spearman correlation for each tissue
  Tissue_correlations <- sapply(Final_expr_list, function(tbl) {
    if (nrow(tbl) < 2){
      NA  # Not enough data to compute correlation
    }
    else{
      cor(tbl$RNA, tbl$protein, method = "spearman", use = "complete.obs")
    }
  })

  Tissue_correlations_tbl <- tibble(tissue = names(Tissue_correlations), Spearman_correlation = Tissue_correlations)

  # Generating a plot using ggplot2 to visualize the Spearman correlation of each tissue
    output_plot <-ggplot(Tissue_correlations_tbl, aes(x = tissue, y = Spearman_correlation, fill = Spearman_correlation)) +
    geom_col() +
    geom_text(aes(label = round(Spearman_correlation, 2)), vjust = -0.5) +
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0,
      name = "Spearman Correlation\n(red = high negative,\nwhite = ~0, blue = high positive)"
    ) +
    labs(
      title = "Spearman Correlation between RNA and Protein Expression \n of the top expressed genes in tissues of interest",
      x = "Tissues",
      y = "Spearman Correlation Values"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )


  # Returning an object output which contains the final plot we created
  output <- list(per_tissue_plot = output_plot)

  return (output)

}







#' Generate Correlation Plots for given list of genes and tissues
#'
#' This function serves as a wrapper to compute and visualize the Spearman correlation
#' between RNA and protein expression between a user's gene list **AND** tissue list of interest.
#' Depending on the user's choice, it calls either `per_gene_plot` or `per_tissue_plot` internally.
#' `per_gene_plot` visualizes the Spearman correlation and plot per gene from the input list.
#' `per_tissue_plot` visualizes the Spearman correlation and plot per tissue from the input list.
#'
#' @note
#' This function computes correlations for a list of genes **AND** a list of tissues.
#' If you wish to input  **ONLY** a list of genes **OR** list of tissues, utilize function `compute_correlation`, or
#' see `help` documentation `?compute_correlation`
#'
#' @param gene_NAMES A character vector of gene symbols. Must contain at least five elements if provided.
#' Valid gene symbols for the \code{gene_NAMES} argument can be accessed from
#' \code{gene_symbols_list}.
#'
#' @param tissue_NAMES A character vector of tissue names. Must contain at least five elements if provided.
#' Valid tissue names as for the \code{tissue_NAMES} argument can be accessed from
#' \code{tissue_map$protein_tissue}
#'
#' @param plot_choice A character string specifying which plot to generate.
#'   Options are `"per_gene"` to visualize correlation for each gene, or `"per_tissue"`
#'   to visualize correlation for each tissue.
#'
#' @return Depending on `plot_choice`, the function returns:
#' \describe{
#'   \item{per_gene_plot}{A plot object showing Spearman correlations per gene (if `plot_choice = "per_gene"`).}
#'   \item{per_tissue_plot}{A plot object showing Spearman correlations per tissue (if `plot_choice = "per_tissue"`).}
#' }
#'
#' @details
#' The function validates the inputs before generating plots:
#' \itemize{
#'   \item If both `gene_NAMES` and `tissue_NAMES` are NULL, it stops execution and throws an error.
#'   \item If both lists have fewer than five elements, it stops execution and throws an error.
#'   \item If only one list is provided, it stops execution and throws an error.
#'   \item If either list has duplicates, it stops execution and throws an error.
#'   \item The `plot_choice` argument must be either `"per_gene"` or `"per_tissue"`. Any other value stops execution and throws an error.
#'   \item Plot key - Spearman correlation values are represented by color intensity. High positive correlations are indicated in dark blue,
#' low correlations around zero are white, and high negative correlations are indicated in dark red.
#' The shade intensity reflects the magnitude of the correlation.
#' }
#'
#' The function then calls:
#' \itemize{
#'   \item `per_gene_plot(tissue_NAMES, gene_NAMES)` if `plot_choice = "per_gene"`.
#'   \item `per_tissue_plot(tissue_NAMES, gene_NAMES)` if `plot_choice = "per_tissue"`.
#' }
#'
#' @note ExprCompareR retrieves data from HPA (Thul et al., 2018) and GTEx (Lonsdale et al., 2013)
#'via packages, HPAanalyze(Tran et al., 2019) and gtexr(Warwick et al., 2025).
#'Because these functions make external API calls, some operations may take time to complete
#'typically around 30 seconds to 1 minute.
#'
#' @examples
#' \dontrun{
#' # Example 1: Generate correlation plot per gene
#'
#' tissues = c("lung", "spleen", "liver", "ovary", "testis")
#' genes = c("MYC", "TP53", "BRCA1", "CRP", "EGFR")
#'
#' result_gene <- correlation_genes_tissues(gene_NAMES = genes, tissue_NAMES = tissues, plot_choice = "per_gene")
#' result_gene$per_gene_plot
#'
#' # Example 2: Generate correlation plot per tissue
#'
#'tissues = c("lung", "spleen", "liver", "ovary", "testis")
#' genes = c("MYC", "TP53", "BRCA1", "CRP", "EGFR")
#'
#' result_gene <- correlation_genes_tissues(gene_NAMES = genes,tissue_NAMES = tissues, plot_choice = "per_tissue")
#' result_gene$per_tissue_plot
#' }
#'
#'@references
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
#'Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). "gtexr: A
#'convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API."
#'Journal of Open Source Software, 10(109), 8249. ISSN 2475-9066,
#'\url{doi:10.21105/joss.08249}, gigs v0.2.1.
#'
#'Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#'Springer-Verlag New York. ISBN 978-3 319-24277-4, \url{https://ggplot2.tidyverse.org.}
#'
#'Wickham H, Fracois R, Henry L, Muller K, Vaughan D (2025).
#'dplyr: A Grammar of Data Manipulation. R package version 1.1.4,
#'\url{https://dplyr.tidyverse.org.}
#'
#'Wickham H, Henry L (2025). purrr: Functional Programming Tools.
#'R package version 1.1.0, \url{https://purrr.tidyverse.org.}
#'
#'Wickham H, Vaughan D, Girlich M (2025). tidyr: Tidy Messy Data.
#'R package version 1.3.1, \url{https://tidyr.tidyverse.org.}
#'
#'
#' @export

correlation_genes_tissues<- function(gene_NAMES, tissue_NAMES, plot_choice){

  # Case 1: Inputs, gene_NAMES and tissue_NAMES are both NULL
  if (is.null(gene_NAMES) && is.null(tissue_NAMES)){
    stop("Please provide at least five gene names AND at least five tissue names")
  }

  # Case 2: Only one list provided (function requires both lists)
  if (is.null(gene_NAMES) || is.null(tissue_NAMES)){
    stop("This function requires both gene_NAMES and tissue_NAMES. For only one list, use compute_correlation().")
  }

  # Case 3: if gene list <5
  if (length(gene_NAMES) < 5){
    stop("gene_NAMES must contain at least five elements.")
  }

  # Case 4: tissue list <5
  if (length(tissue_NAMES) < 5){
    stop("tissue_NAMES must contain at least five elements.")
  }

  # Case 5: Check for duplicates in gene_NAMES
  if (any(duplicated(gene_NAMES))) {
    dupes <- unique(gene_NAMES[duplicated(gene_NAMES)])
    stop(
      paste0(
        "Duplicate gene names detected: ",
        paste(dupes, collapse = ", "),
        ".\nPlease remove duplicates before proceeding."
      )
    )
  }

  # Case 6: Check for duplicates in tissue_NAMES
  if (any(duplicated(tissue_NAMES))) {
    dupes <- unique(tissue_NAMES[duplicated(tissue_NAMES)])
    stop(
      paste0(
        "Duplicate tissue names detected: ",
        paste(dupes, collapse = ", "),
        ".\nPlease remove duplicates before proceeding."
      )
    )
  }


  # Validate plot_choice
  if (!plot_choice %in% c("per_gene","per_tissue")){
    stop("Choice of variable plot_choice is invalid. Please choose
            between per_gene and per_tissue")
  }

  # If plot_choice is per_gene, call function per_gene_plot
  if (plot_choice == "per_gene"){
    output<- per_gene_plot(gene_NAMES, tissue_NAMES)
    return(output)
  }

  # If plot_choice is per_tissue, call function per_tissue_plot
  else if (plot_choice == "per_tissue"){
    output<- per_tissue_plot(gene_NAMES, tissue_NAMES)
    return(output)
  }


}




#' Generate Spearman Correlation Plot Per Gene Across User's Gene List and Tissue List Of Interest
#'
#' This function calculates the Spearman correlation between RNA and protein expression
#' for each gene in \code{gene_NAMES} across the specified tissues in \code{tissue_NAMES}.
#' It produces a bar plot visualizing the correlation for each gene.
#'
#' @param gene_NAMES A character vector of gene symbols. Must contain at least five elements if provided.
#' Valid gene symbols for the \code{gene_NAMES} argument can be accessed from
#' \code{gene_symbols_list}.
#'
#' @param tissue_NAMES A character vector of tissue names. Must contain at least five elements if provided.
#' Valid tissue names as for the \code{tissue_NAMES} argument can be accessed from
#' \code{tissue_map$protein_tissue}

#'
#' @return A list containing:
#' \describe{
#'   \item{per_gene_plot}{A \code{ggplot2} object showing Spearman correlation values
#'   for each gene across the specified tissues.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Converts tissue names to GTEx-compatible names using \code{convert_to_gtex}.
#'   \item Converts gene symbols to GENCODE IDs using \code{get_gtex_gencode_ids}.
#'   \item Queries RNA expression data from GTEx for each gene and tissue combination.
#'   \item Aggregates RNA expression by averaging multiple samples per gene per tissue.
#'   \item Downloads normal tissue protein expression data using \code{HPAanalyze::hpaDownload}.
#'   \item Filters protein data to include only the genes and tissues of interest.
#'   \item Handles missing gene-tissue combinations, filling missing entries with \code{NA} and issuing a warning.
#'   \item Converts protein expression levels to numeric using \code{protein_expr_values}.
#'   \item Aggregates protein data by mean per gene and tissue.
#'   \item Computes Spearman correlations between RNA and protein expression for each gene.
#'   \item Generates a bar plot showing correlations for all genes.
#'
#'
#'  - Plot key - Spearman correlation values are represented by color intensity. High positive correlations are indicated in dark blue,
#' low correlations around zero are white, and high negative correlations are indicated in dark red.
#' The shade intensity reflects the magnitude of the correlation.
#' }
#'
#'
#'@note ExprCompareR retrieves data from HPA (Thul et al., 2018) and GTEx (Lonsdale et al., 2013)
#'via packages, HPAanalyze(Tran et al., 2019) and gtexr(Warwick et al., 2025).
#'Because these functions make external API calls, some operations may take time to complete
#'typically around 30 seconds to 1 minute.
#'
#' @examples
#' \dontrun{
#' # Example 1: Generate correlation plot per gene
#'
# tissues= c("lung", "spleen", "liver", "ovary", "testis")
# genes = c("MYC", "TP53", "BRCA1", "CRP", "EGFR")
#'
#' result_gene <- per_gene_plot(gene_NAMES = genes, tissue_NAMES = tissues)
#' result_gene$per_gene_plot
#'}
#'
#' @references
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
#'Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). "gtexr: A
#'convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API."
#'Journal of Open Source Software, 10(109), 8249. ISSN 2475-9066,
#'\url{doi:10.21105/joss.08249}, gigs v0.2.1.
#'
#'Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#'Springer-Verlag New York. ISBN 978-3 319-24277-4, \url{https://ggplot2.tidyverse.org.}
#'
#'Wickham H, Francois R, Henry L, Muller K, Vaughan D (2025).
#'dplyr: A Grammar of Data Manipulation. R package version 1.1.4,
#'\url{https://dplyr.tidyverse.org.}
#'
#'Wickham H, Henry L (2025). purrr: Functional Programming Tools.
#'R package version 1.1.0, \url{https://purrr.tidyverse.org/.}
#'
#'Wickham H, Vaughan D, Girlich M (2025). tidyr: Tidy Messy Data.
#'R package version 1.3.1, \url{https://tidyr.tidyverse.org.}
#'
#'
#'
#' @import HPAanalyze
#' @import gtexr
#' @import ggplot2
#' @import tidyr
#' @import purrr
#' @import dplyr
#' @importFrom stats cor median
#' @importFrom utils head

#' @keywords internal

per_gene_plot <- function(gene_NAMES, tissue_NAMES){

  # Query RNA expression

  # Convert the tissue_NAMES to RNA compatible names using our helper function-
  # convert_to_gtex
  tissue_list_RNA <- convert_to_gtex(tissue_NAMES)
  # Convert gene symbols to RNA compatible GENCODE IDs using our helper function-
  # get_gtex_gencode_ids
  gene_ids <- get_gtex_gencode_ids(gene_NAMES)


  # Querying gtexr to obtain RNA expression for the given gene and tissue list
  RNA_expr<- suppressMessages(gtexr::get_gene_expression(gencodeIds = gene_ids, tissueSiteDetailIds = tissue_list_RNA))


  # Computing the average of different samples of a gene in the tissue
  # as we have multiple samples of a given gene.
  RNA_expr <- RNA_expr %>%
    mutate(avg_expr = map_dbl(data, ~ mean(.x, na.rm = TRUE)))

  # Query Protein expression

  # # Using HPAanalyze's built-in function hpaDownload to get protein expression data
  # downloadedData <- HPAanalyze::hpaDownload(downloadList='histology', version='example')
  # # Subsetting downloadedData to get normal tissue data
  # normal_tissue <- downloadedData$normal_tissue


  # Filter protein data to only include genes from gene_NAMES and tissues in
  # tissue_NAMES
  protein_levels <- normal_tissue %>%
    dplyr::filter(gene %in% gene_NAMES, tissue %in% tissue_NAMES) %>%
    dplyr::select(gene, tissue, level)

  # Check if there is a gene not found in a certain tissue

  # Summarize counts per gene-tissue and suppress output
  invisible(protein_levels %>%
              dplyr::group_by(gene, tissue) %>%
              dplyr::summarise(n = n(), .groups = "drop"))

  # Create all possible gene-tissue combinations
  all_combos <- expand.grid(gene = gene_NAMES, tissue = tissue_NAMES, stringsAsFactors = FALSE)

  # Identify missing combinations
  missing_combos <- dplyr::anti_join(all_combos, protein_levels, by = c("gene", "tissue"))

  # Check for missing combos and issue a warning if any exist
  if (nrow(missing_combos) != 0) {
    warning(
      paste(
        "Some gene-tissue combinations are missing in protein data. Check 'missing_combos' tibble for details.",
        "These missing combos indicate that specific genes of interest are not present in the",
        "tissues of interest.",
        "You can either remove these genes from your analysis, replace these values or they will be filled by NA values.",
        sep = "\n"
      )
    )

    missing_combos <- missing_combos %>%
      dplyr::mutate(level = NA)


    # Append missing rows to protein_levels
    protein_levels <- dplyr::bind_rows(protein_levels, missing_combos)
  }


  # Converting the protein expression string values to numeric values
  protein_levels<- protein_levels %>%
    mutate(level = protein_expr_values(level))


  # Rename columns for consistency
  RNA_expr_clean <- RNA_expr %>%
    rename(tissue = tissueSiteDetailId, gene = geneSymbol, RNA = avg_expr)

  protein_levels_clean <- protein_levels %>%
    rename(protein = level)

  # Aggregate protein by gene and tissue (mean)
  protein_summary <- protein_levels_clean %>%
    dplyr::group_by(gene, tissue) %>%
    dplyr::summarise(protein = mean(protein, na.rm = TRUE), .groups = "drop")


  # Convert RNA tissues to protein
  RNA_expr_clean <- RNA_expr_clean %>%
    dplyr::left_join(tissue_map, by = c("tissue" = "RNA_tissue")) %>%
    dplyr::select(-tissue) %>%
    dplyr::rename(tissue = protein_tissue)


  # Create a tibble for each gene, in which every row is the tissue from tissue_NAMES.
  # Each tibble has two columns, RNA and protein, and each value represents that gene's
  # RNA and protein expression in that tissue.
  gene_tables <- purrr::map(unique(RNA_expr_clean$gene), function(g) {

    RNA_sub <- RNA_expr_clean %>%
      dplyr::filter(gene == g) %>%
      dplyr::select(tissue, RNA)

    protein_sub <- protein_summary %>%
      dplyr::filter(gene == g) %>%
      dplyr::select(tissue, protein)

    joined <- dplyr::inner_join(RNA_sub, protein_sub, by = "tissue")

    return(joined)
  })

  names(gene_tables) <- unique(RNA_expr_clean$gene)

  # Calculate Spearman correlation for each gene table
  gene_correlations <- purrr::map_dbl(gene_tables, function(joined) {
    cor(joined$RNA, joined$protein, use = "complete.obs", method = "spearman")
  })

  # Convert to a tibble for easy viewing
  gene_correlations_tbl <- tibble::tibble(
    gene = names(gene_correlations),
    spear_corr = gene_correlations
  )


  # Generating a plot using ggplot2 to visualize the Spearman correlation of each gene.

  output_plot <- ggplot(gene_correlations_tbl, aes(x = gene, y = spear_corr, fill = spear_corr)) +
    geom_col() +
    geom_text(aes(label = round(spear_corr, 2), vjust = -0.5)) +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0,  name = "Spearman Correlation\n(red = high negative,\nwhite = ~0, blue = high positive)") +
    labs(
      title = "Spearman Correlation between RNA and Protein Expression \n of genes of interest across tissues of interest",
      x = "Genes",
      y = "Spearman Correlation Value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )

  # Returning an object output which contains the final plot we created
  output<- list(per_gene_plot = output_plot)

  return(output)

}




#' Generate Spearman Correlation Plot Per Tissue Across User's Gene List and Tissue List Of Interest
#'
#' This function calculates the Spearman correlation between RNA and protein expression
#' for the top genes in \code{gene_NAMES} across the specified tissues in \code{tissue_NAMES}.
#' It produces a bar plot visualizing the correlation for each tissue.
#'
#'@param gene_NAMES A character vector of gene symbols. Must contain at least five elements if provided.
#' Valid gene symbols for the \code{gene_NAMES} argument can be accessed from
#' \code{gene_symbols_list}.
#'
#' @param tissue_NAMES A character vector of tissue names. Must contain at least five elements if provided.
#' Valid tissue names as for the \code{tissue_NAMES} argument can be accessed from
#' \code{tissue_map$protein_tissue}
#'
#' @return A list containing:
#' \describe{
#'   \item{per_tissue_plot}{A \code{ggplot2} object showing Spearman correlation values
#'   for each tissue across the specified genes.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Converts tissue names to GTEx-compatible names using \code{convert_to_gtex}.
#'   \item Converts gene symbols to GENCODE IDs using \code{get_gtex_gencode_ids}.
#'   \item Queries RNA expression data from GTEx for each gene and tissue combination.
#'   \item Aggregates RNA expression by averaging multiple samples per gene per tissue.
#'   \item Downloads normal tissue protein expression data using \code{HPAanalyze::hpaDownload}.
#'   \item Filters protein data to include only the genes and tissues of interest.
#'   \item Handles missing gene-tissue combinations by filling missing entries with \code{NA} and issuing a warning.
#'   \item Converts protein expression levels to numeric using \code{protein_expr_values}.
#'   \item Aggregates protein data by mean per gene per tissue.
#'   \item Computes Spearman correlations between RNA and protein expression for each tissue.
#'   \item Generates a bar plot showing correlations for all tissues.
#'
#'  - Plot key - Spearman correlation values are represented by color intensity. High positive correlations are indicated in dark blue,
#' low correlations around zero are white, and high negative correlations are indicated in dark red.
#' The shade intensity reflects the magnitude of the correlation.
#' }
#'
#'@note ExprCompareR retrieves data from HPA (Thul et al., 2018) and GTEx (Lonsdale et al., 2013)
#'via packages, HPAanalyze(Tran et al., 2019) and gtexr(Warwick et al., 2025).
#'Because these functions make external API calls, some operations may take time to complete
#'typically around 30 seconds to 1 minute.
#'
#' @examples
#' \dontrun{
#' # Example 2: Generate correlation plot per tissue
#'
# tissues= c("lung", "spleen", "liver", "ovary", "testis")
# genes = c("MYC", "TP53", "BRCA1", "CRP", "EGFR")
#'
#' result_gene <- per_tissue_plot(gene_NAMES = genes,tissue_NAMES = tissues)
#' result_gene$per_tissue_plot
#' }
#'
#' @references
#'
#'Müller K, Wickham H (2025). tibble: Simple Data Frames.
#' R package version 3.3.0, \url{https://tibble.tidyverse.org/.}
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
#'Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). "gtexr: A
#'convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API.
#'Journal of Open Source Software, 10(109), 8249. ISSN 2475-9066,
#'\url{doi:10.21105/joss.08249}, gigs v0.2.1.
#'
#'Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#'Springer-Verlag New York. ISBN 978-3 319-24277-4,
#'\url{https://ggplot2.tidyverse.org.}
#'
#'Wickham H, Francois R, Henry L, Muller K, Vaughan D (2025).
#'dplyr: A Grammar of Data Manipulation. R package version 1.1.4,
#'\url{https://dplyr.tidyverse.org.}
#'
#'Wickham H, Henry L (2025). purrr: Functional Programming Tools.
#'R package version 1.1.0, \url{https://purrr.tidyverse.org/.}
#'
#'Wickham H, Vaughan D, Girlich M (2025). tidyr: Tidy Messy Data.
#'R package version 1.3.1, \url{https://tidyr.tidyverse.org.}
#'
#'
#' @import HPAanalyze
#' @import gtexr
#' @import ggplot2
#' @import tidyr
#' @import purrr
#' @import dplyr
#' @import tibble
#' @importFrom stats cor median
#' @importFrom utils head

#' @keywords internal

per_tissue_plot <- function (gene_NAMES, tissue_NAMES){

  # Query RNA expression

  # Convert the tissue_NAMES to RNA compatible names using our helper function-
  # convert_to_gtex
  tissue_list_RNA <- convert_to_gtex(tissue_NAMES)
  # Convert gene symbols to RNA compatible GENCODE IDs using our helper function-
  # get_gtex_gencode_ids
  gene_ids <- get_gtex_gencode_ids(gene_NAMES)

  # Querying gtexr to obtain RNA expression for the given gene and tissue list
  RNA_expr <- suppressMessages(gtexr::get_gene_expression(gencodeIds = gene_ids, tissueSiteDetailIds = tissue_list_RNA))
  head(RNA_expr)

  # Computing the average of different samples of a gene in the tissue
  # as we have multiple samples of a given gene.
  RNA_expr <- RNA_expr %>%
    mutate(avg_expr = map_dbl(data, ~ mean(.x, na.rm = TRUE)))


  # Query Protein Data

  # Filter protein data to only include genes from gene_NAMES and tissues in
  # tissue_NAMES
  protein_levels <- normal_tissue %>%
    dplyr::filter(gene %in% gene_NAMES, tissue %in% tissue_NAMES) %>%
    dplyr::select(gene, tissue, level)


  # Summarise counts per gene-tissue and suppress output
  invisible(protein_levels %>%
              dplyr::group_by(gene, tissue) %>%
              dplyr::summarise(n = n(), .groups = "drop"))

  # Check if there is a gene not found in a certain tissue

  # Create all possible gene-tissue combinations
  all_combos <- expand.grid(gene = gene_NAMES, tissue = tissue_NAMES, stringsAsFactors = FALSE)

  # Identify missing combinations
  missing_combos <- dplyr::anti_join(all_combos, protein_levels, by = c("gene", "tissue"))

  # Check for missing combos and issue a warning if they exist
  if (nrow(missing_combos) != 0) {
    warning(
      paste(
        "Some gene-tissue combinations are missing in protein data. Check 'missing_combos' tibble for details.",
        "These missing combos indicate that specific genes of interest are not present in the",
        "tissues of interest.",
        "You can either remove these genes from your analysis, replace these values or they will be filled by NA values.",
        sep = "\n"
      )
    )

    missing_combos <- missing_combos %>%
      dplyr::mutate(level = NA)


    # Append missing rows to protein_levels
    protein_levels <- dplyr::bind_rows(protein_levels, missing_combos)
  }


  # Converting the protein expression string values to numeric values
  protein_levels<- protein_levels %>%
    mutate(level = protein_expr_values(level))


  # Rename columns for consistency
  RNA_expr_clean <- RNA_expr %>%
    rename(tissue = tissueSiteDetailId, gene = geneSymbol, RNA = avg_expr)

  protein_levels_clean <- protein_levels %>%
    rename(protein = level)

  # Aggregate protein by gene and tissue (mean)
  protein_summary <- protein_levels_clean %>%
    dplyr::group_by(gene, tissue) %>%
    dplyr::summarise(protein = mean(protein, na.rm = TRUE), .groups = "drop")


  # Convert RNA tissues to protein
  RNA_expr_clean <- RNA_expr_clean %>%
    dplyr::left_join(tissue_map, by = c("tissue" = "RNA_tissue")) %>%
    dplyr::select(-tissue) %>%
    dplyr::rename(tissue = protein_tissue)


  # Generating a tibble for each tissue, with each row being the gene's RNA
  # and protein expression values
  Final_tibble <- lapply(unique(RNA_expr_clean$tissue), function(tissue_name) {

    # Clean RNA table
    RNA_tbl <- RNA_expr_clean %>%
      filter(tissue == tissue_name) %>%
      mutate(
        gene = trimws(as.character(gene)),
        RNA = as.numeric(RNA)
      ) %>%
      select(gene, RNA)

    # Clean protein table
    protein_tbl <- protein_summary %>%
      mutate(
        gene = trimws(as.character(gene)),
        tissue = trimws(as.character(tissue)),
        protein = as.numeric(protein)
      ) %>%
      filter(tissue == tolower(tissue_name)) %>%
      select(gene, protein)

    # Merge RNA and protein
    merged_tbl <- left_join(RNA_tbl, protein_tbl, by = "gene") %>%
      mutate(protein = ifelse(is.na(protein), 0, protein))

    return(merged_tbl)
  })

  names(Final_tibble) <- unique(RNA_expr_clean$tissue)


  # Calculate Spearman Coefficient for each tissue
  Tissue_correlations <- sapply(Final_tibble, function(tbl) {
    if (nrow(tbl) < 2){
      NA  # Not enough data to compute correlation
    }
    else{
      cor(tbl$RNA, tbl$protein, method = "spearman", use = "complete.obs")
    }
  })

  Tissue_correlations_tbl <- tibble(tissue = names(Tissue_correlations), Spearman_correlation = Tissue_correlations)


  # Generating a plot using ggplot2 to visualize the Spearman correlation of each tissue.

  output_plot <- ggplot(Tissue_correlations_tbl, aes(x = tissue, y = Spearman_correlation, fill = Spearman_correlation)) +
    geom_col() +
    geom_text(aes(label = round(Spearman_correlation, 2)), vjust = -0.5) +
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0,
      name = "Spearman Correlation\n(red = high negative,\nwhite = ~0, blue = high positive)"
    )+
    labs(
      title = "Spearman Correlation between RNA and Protein Expression \n of genes of interest across tissues of interest",
      x = "Tissues",
      y = "Spearman Correlation Values "
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )

  # Returning an object output which contains the final plot we created
  output <- list(per_tissue_plot = output_plot)

  return (output)

}


# [END]
