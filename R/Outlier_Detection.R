utils::globalVariables(c(
  "tissue", "reliability", "level", "ensembl", "gene",
  "rank_value", "normal_rank", "data", "geneSymbol",
  "avg_expr", "RNA", "Protein", "RNA_log", "is_outlier"
))


#'Detects outlier genes present in a given tissue of interest, which have a
#'large delta between their RNA expression and protein expression values
#'
#'The function given a tissue of interest, compares the RNA expression and protein expression
#'of all genes present in the given tissue, and identifies genes which have a large difference
#'in their RNA and protein expression. These genes are identified as outlier genes, and are most
#'likely candidates for post transcriptional modifications. This function returns the list of outlier genes
#'as well as a plot highlighting the outlier genes.
#'
#'@param input_tissue Character String specifying the tissue of interest to analyze.
#'The name must be present in the given list of tissues:
#' Valid tissue names are listed in \code{names(tissue_map)}.
#' Any name not in this list will issue a warning.
#'
#'
#'@return Returns an object of class list with two components:
#'\itemize{
#'\item \code{outlier_vector}: A character vector of outlier genes names detected
#'\item \code{outlier_plot}: A ggplot object highlighting the outlier genes in red
#'}
#'
#'@details
#'The function downloads protein expression data from the Human Protein Atlas
#' using \pkg{HPAanalyze} and RNA expression data from the GTEx portal using \pkg{gtexr}.
#' It compares both modalities for the specified tissue and flags outlier genes
#' based on interquartile range (IQR) thresholds for RNA and protein expression.
#'
#'@examples
#'\dontrun{
#' # Example 1: Detect outliers in lung tissue
#'
#' output <- detect_outliers(input_tissue = "lung")
#'
#' output$outlier_vector
#' output$outlier_plot
#'
#'
#' # Example 2: Detect outlier in the ovary tissue
#'
#'output <- detect_outliers(input_tissue = "ovary")
#'
#'output$outlier_vector
#'output$outlier_plot
#'}
#'
#'@references
#'
#'Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). “gtexr: A
#'convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API.”
#'Journal of Open Source So ware, 10(109), 8249. ISSN
#'2475-9066, doi:10.21105/joss.08249, gigs v0.2.1.
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'“HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data.” BMC
#'
#'Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#'Springer-Verlag New York. ISBN 978-3 319-24277-4, https://ggplot2.tidyverse.org.
#'
#'Wickham H, Henry L (2025). purrr: Functional Programming Tools.
#'R package version 1.1.0, https://purrr.tidyverse.org/.
#'
#'Wickham H, Vaughan D, Girlich M (2025). tidyr: Tidy Messy Data.
#'R package version 1.3.1, https://tidyr.tidyverse.org.
#'
#'Wickham H, François R, Henry L, Müller K, Vaughan D (2025).
#'dplyr: A Grammar of Data Manipulation. R package version 1.1.4,
#'https://dplyr.tidyverse.org.
#'
#'
#'@import HPAanalyze
#'@import gtexr
#'@import ggplot2
#'@import dplyr
#'@import purrr
#'@import tidyr
#'@importFrom stats quantile
#'@importFrom utils data

#'@export
detect_outliers <- function(input_tissue){

  # Using HPAanalyze's built-in function hpaDownload to get protein expression data
  downloadedData <- hpaDownload(downloadList='histology', version='example')
  # Subsetting downloadedData to get normal tissue data
  normal_tissue <- downloadedData$normal_tissue


  # Subset normal_tissue to obtain protein expression of our tissue of interest,
  # and only consider values which have proper reliability to filter out protein
  # expression values that have low reliability
  tissue_data <- normal_tissue %>%
    filter(tissue == input_tissue,
           reliability %in% c("Enhanced", "Supported"))


  # Converting the protein expression string values to numeric values
  tissue_data <- tissue_data %>%
    dplyr::mutate(rank_value = protein_expr_values(level))


  # Computing the average of different samples of a gene in the tissue
  # as we have multiple samples of a given gene.
  tissue_data <- tissue_data %>%
    dplyr::group_by(ensembl, gene) %>%
    dplyr::summarise(normal_rank = mean(rank_value, na.rm = TRUE), .groups = "drop")


  # Filtering to only keep genes who have a expression > 0 to eliminate gene's which
  # have 0.0 protein expression
  tissue_data <- tissue_data %>% filter(normal_rank > 0)

  # Extracting the names of the genes we filtered in protein expression
  gene_names <- unique(tissue_data$gene)


  # To query for RNA expression of these genes, we convert these gene names to
  # GENCODE IDs as the package gtexr takes GENCODE IDs as input.

  # As there are many genes, we query them in batches of 500
  # to avoid URL fetching failure
  batch_size <- 500
  all_ids <- c()

  for (i in seq(1, length(gene_names), by = batch_size)) {
    batch <- gene_names[i:min(i + batch_size - 1, length(gene_names))]
    batch_ids <- get_gtex_gencode_ids(batch)
    all_ids <- c(all_ids, batch_ids)
    Sys.sleep(0.5)  # pause to avoid API rate limiting
  }

  # Filter out those genes whose GENCODE IDs could not be found.
  all_ids <- all_ids[!is.na(all_ids)]


  # Convert the input_tissue name to the correct tissue name gtexr uses.
  RNA_input_tissue = convert_to_rna(input_tissue)


  # Query RNA expression for those genes found in the input tissue
  # Split IDs into batches of 100 to avoid URL fetching failure
  batch_size <- 100
  batches <- split(all_ids, ceiling(seq_along(all_ids)/batch_size))

  RNA_result_list <- list()

  for (i in seq_along(batches)) {
    cat("Fetching batch", i, "of", length(batches), "\n")
    RNA_result_list[[i]] <- gtexr::get_gene_expression(
      gencodeIds = batches[[i]],
      tissueSiteDetailIds = RNA_input_tissue
    )
  }

  # Combine results into one dataframe
  RNA_result <- do.call(rbind, RNA_result_list)


  # Average out the RNA expression of various samples of a given gene, to obtain
  # obtain one single RNA expression value
  RNA_result <- RNA_result %>%
    mutate(avg_expr = map_dbl(data, ~ mean(.x, na.rm = TRUE)))


  # Create a table in which each row represents a gene, and its respective RNA
  # and protein expression
  expr_table <- data.frame(
    gene = gene_names
  ) %>%
    left_join(RNA_result %>% select(geneSymbol, avg_expr), by = c("gene" = "geneSymbol")) %>%
    left_join(tissue_data %>% select(gene, normal_rank), by = "gene") %>%
    rename(RNA = avg_expr, Protein = normal_rank)


  # Removing NA values from the table
  expr_table <- expr_table %>%
    tidyr::drop_na(RNA, Protein)


  # Find the outliers from the genes found in the tissue of interest.

  #  Log-transform only RNA values, as they are far right skewed as compared to
  # the protein values to have a fair comparison
  expr_table <- expr_table %>%
    mutate(RNA_log = log2(RNA + 1))  # +1 to avoid log(0)

  # Compute IQR thresholds to identify outlier genes
  rna_Q1 <- quantile(expr_table$RNA_log, 0.25)
  rna_Q3 <- quantile(expr_table$RNA_log, 0.75)
  rna_IQR <- rna_Q3 - rna_Q1
  rna_lower <- rna_Q1 - 1.5 * rna_IQR
  rna_upper <- rna_Q3 + 1.5 * rna_IQR

  protein_Q1 <- quantile(expr_table$Protein, 0.25)
  protein_Q3 <- quantile(expr_table$Protein, 0.75)
  protein_IQR <- protein_Q3 - protein_Q1
  protein_lower <- protein_Q1 - 1.5 * protein_IQR
  protein_upper <- protein_Q3 + 1.5 * protein_IQR

  # Identify outliers between RNA log vs Protein raw values
  expr_table <- expr_table %>%
    mutate(is_outlier = (RNA_log < rna_lower | RNA_log > rna_upper |
                           Protein < protein_lower | Protein > protein_upper))

  outlier_genes <- expr_table %>%
    filter(is_outlier) %>%
    pull(gene)


  # Plot the RNA expression vs Protein expression for all these genes.
  # Specifically highlight the genes which are outliers.
  plots <- ggplot(expr_table, aes(x = RNA_log, y = Protein)) +
    geom_point(aes(color = is_outlier), alpha = 0.7) +
    geom_text(aes(label = ifelse(is_outlier, gene, "")),
              vjust = -0.5, size = 3) +
    scale_color_manual(values = c("black", "red")) +
    labs(x = "Log2 RNA expression", y = "Protein expression",
         title = paste("RNA vs Protein Correlation in ", RNA_input_tissue, "\n Outlier genes indicated in red"),
         color = "Outlier") +
    theme_minimal()


  # return an object which has the outlier genes list, and the plot
  output <- list(outlier_vector = outlier_genes, outlier_plot = plots)

  return (output)
}

# [END]
