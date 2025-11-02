#' Tissue mapping: HPAanalyze to gtexr
#'
#' A named character vector that maps tissue names from the Human Protein Atlas (HPAanalyze)
#' to the corresponding tissue names used in the GTEx portal via the gtexr package.
#' This is useful for converting HPA tissue names to GTEx-compatible names
#' when querying RNA expression data.
#'
#'@format A named character vector of length 37.
#' \describe{
#'   \item{names}{HPAanalyze tissue names (character)}
#'   \item{values}{Corresponding GTEx tissue names (character)}
#' }
#'
#'@source  Internal dataset created for ExprCompareR package.
#'Data from Human Protein Atlas (\url{https://www.proteinatlas.org/}) via HPAanalyze
#'  and GTEx Portal (\url{https://gtexportal.org/}) via gtexr
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
#'GTEx Consortium (2013). The Genotype-Tissue Expression (GTEx) project.
#'Nature genetics, 45(6), 580–585. https://doi.org/10.1038/ng.2653
#'
#'Thul PJ, Lindskog C. The human protein atlas: A spatial map of the human proteome.
#'Protein Sci. 2018 Jan;27(1):233-244. doi: 10.1002/pro.3307.
#'Epub 2017 Oct 10. PMID: 28940711; PMCID: PMC5734309.
#'
#'@examples
#' data(tissue_map)
#' length(tissue_map)
#' head(tissue_map)
#' names(tissue_map)
#' unname(tissue_map)
#' \dontrun{
#' tissue_map
#' }
"tissue_map"
