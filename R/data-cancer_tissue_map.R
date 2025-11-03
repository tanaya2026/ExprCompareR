#'Mapping of Cancer Types to Corresponding Normal Tissues
#'
#'This dataset provides a mapping between various cancer types and the normal
#' tissues from which these cancers originate or are typically associated. It can
#' be used to translate cancer type labels into corresponding tissue labels for
#' downstream analyses, such as comparing RNA and protein expression between
#' cancers and their normal tissues.
#'
#' @format A tibble with 20 rows and 2 columns:
#' \describe{
#'   \item{cancer}{Character string. The name of the cancer type.}
#'   \item{tissue}{Character string. The corresponding normal tissue.}
#' }
#'@source Internal dataset created for ExprCompareR package.
#'Data from Human Protein Atlas (\url{https://www.proteinatlas.org/}) via HPAanalyze
#'
#' @references
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
#'@examples
#'data(cancer_tissue_map)
#'length(cancer_tissue_map)
#' head(cancer_tissue_map)
#' cancer_tissue_map$cancer
#' cancer_tissue_map$tissue
#' \dontrun{
#' cancer_tissue_map
#' }
"cancer_tissue_map"
