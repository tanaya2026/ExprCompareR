#' Normal tissue protein expression data
#'
#' A data frame containing protein expression levels across various human tissues
#' and cell types, retrieved from the Human Protein Atlas.
#'
#' @format A data frame with 6 variables:
#' \describe{
#'   \item{ensembl}{Ensembl gene identifier (character).}
#'   \item{gene}{HGNC-approved gene symbol (character).}
#'   \item{tissue}{Name of the human tissue (character).}
#'   \item{cell_type}{Specific cell type within the tissue (character).}
#'   \item{level}{Protein expression level, one of: `"High"`, `"Medium"`, `"Low"`, or `"Not detected"`.}
#'   \item{reliability}{Confidence level in expression annotation, e.g. `"Enhanced"`, `"Supported"`, or `"Uncertain"`.}
#' }
#'
#' @details
#' This dataset can be used to analyze baseline protein expression in normal
#' human tissues. It is part of the Human Protein Atlas histology data.
#' The dataset `normal_tissue` was derived from the
#' Human Protein Atlas `histology` data using HPAanalyze, which can be downloaded using:
#' \preformatted{
#' downloadedData <- hpaDownload(downloadList = 'histology', version = 'example')
#' normal_tissue <- downloadedData$normal_tissue
#' }
#'
#' @source  Data derived from Human Protein Atlas (\url{https://www.proteinatlas.org/}),
#'via HPAanalyze
#'
#'@references
#'
#'Thul PJ, Lindskog C. The human protein atlas: A spatial map of the human proteome.
#'Protein Sci. 2018 Jan;27(1):233-244. doi: 10.1002/pro.3307.
#'Epub 2017 Oct 10. PMID: 28940711; PMCID: PMC5734309.
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'“HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data.” MC Bioinformatics 20, 463 (2019).
#'https://doi.org/10.1186/s12859-019-3059-z
#'
#' @examples
#' data(normal_tissue)
#' head(normal_tissue)
#' names(normal_tissue)
#' \dontrun{
#' normal_tissue
#' }
"normal_tissue"
