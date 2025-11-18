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
#'
#' @references
#'
#'Thul PJ, Lindskog C. The human protein atlas: A spatial map of the human proteome.
#'Protein Sci. 2018 Jan;27(1):233-244. \url{doi: 10.1002/pro.3307}.
#'Epub 2017 Oct 10. PMID: 28940711; PMCID: PMC5734309.
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
#'
"cancer_tissue_map"




#' Valid list of gene symbols analysis in ExprCompareR
#'
#' A character vector of unique valid gene symbols from the HPA normal tissue dataset
#' used by ExprCompareR.
#' The symbols can be converted to their respective GENCODE IDs for retrieval
#' of their RNA expression data.
#'
#' @format A character vector of length 15313 (one element per gene).
#'
#' @source Derived from Human Protein Atlas (\url{https://www.proteinatlas.org/}),
#' via HPAanalyze
#'
#'
#' @examples
#' data(gene_symbols_list)
#' length(gene_symbols_list)
#' head(gene_symbols_list)
#' \dontrun{
#' gene_symbols_list
#' }
#'
#' @references
#'
#'Thul PJ, Lindskog C. The human protein atlas: A spatial map of the human proteome.
#'Protein Sci. 2018 Jan;27(1):233-244. \url{doi: 10.1002/pro.3307.}
#'Epub 2017 Oct 10. PMID: 28940711; PMCID: PMC5734309.
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
#'
#'
"gene_symbols_list"



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
#'
#' @examples
#' data(normal_tissue)
#' head(normal_tissue)
#' names(normal_tissue)
#' \dontrun{
#' normal_tissue
#' }
#'
#'@references
#'
#'Thul PJ, Lindskog C. The human protein atlas: A spatial map of the human proteome.
#'Protein Sci. 2018 Jan;27(1):233-244. doi: 10.1002/pro.3307.
#'Epub 2017 Oct 10. PMID: 28940711; PMCID: PMC5734309.\url{https://doi.org/10.1002/pro.3307}
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
#'
#'
"normal_tissue"



#' Cancer pathology protein expression data
#'
#' A data frame describing protein expression profiles and prognostic
#' associations across multiple cancer types, as reported in the Human Protein Atlas.
#'
#' @format A data frame with 11 variables:
#' \describe{
#'   \item{ensembl}{Ensembl gene identifier (character).}
#'   \item{gene}{HGNC-approved gene symbol (character).}
#'   \item{cancer}{Name of the cancer type (character).}
#'   \item{high}{Number or proportion of samples with high protein expression.}
#'   \item{medium}{Number or proportion of samples with medium expression.}
#'   \item{low}{Number or proportion of samples with low expression.}
#'   \item{not_detected}{Number or proportion of samples where protein was not detected.}
#'   \item{prognostic_favorable}{Logical or categorical indicator of favorable prognostic association.}
#'   \item{unprognostic_favorable}{Indicator for unprognostic favorable cases.}
#'   \item{prognostic_unfavorable}{Indicator of unfavorable prognostic association.}
#'   \item{unprognostic_unfavorable}{Indicator for unprognostic unfavorable cases.}
#' }
#'
#' @details
#' This dataset can be used to explore relationships between protein expression
#' and cancer prognosis across different tumor types.
#' The dataset `pathology` was derived from the
#' Human Protein Atlas `histology` data, which can be downloaded using:
#' \preformatted{
#' downloadedData <- hpaDownload(downloadList = 'histology', version = 'example')
#' pathology <- downloadedData$pathology
#' }
#'
#' @seealso [HPAanalyze::hpaDownload]
#'
#' @source Data derived from Human Protein Atlas \url{https://www.proteinatlas.org/}
#'
#' @examples
#' data(pathology)
#' head(pathology)
#' names(pathology)
#' \dontrun{
#' pathology
#' }
#'
#' @references
#'
#'Thul PJ, Lindskog C. The human protein atlas: A spatial map of the human proteome.
#'Protein Sci. 2018 Jan;27(1):233-244. doi: 10.1002/pro.3307.
#'Epub 2017 Oct 10. PMID: 28940711; PMCID: PMC5734309. \url{https://doi.org/10.1002/pro.3307}
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
"pathology"



#' Default list of protein tissues used for correlation analysis
#'
#'A character vector of the default tissues that function correlation_per_gene
#'uses to analyze RNA vs protein expression when the input is a list of genes.
#'These tissue names follow the format of HPAanalyze and can only be used to query
# protein expression.
#'
#'@format A character vector of length 10 (one element per tissue)
#'
#'@source Internal dataset created for ExprCompareR package.
#'Data derived from Human Protein Atlas (\url{https://www.proteinatlas.org/}),
#'via HPAanalyze
#'
#'@examples
#'data(tissue_list_protein)
#'length(tissue_list_protein)
#'head(tissue_list_protein)
#'\dontrun{
#'tissue_list_protein
#'}
#'
#'@references
#'
#'Thul PJ, Lindskog C. The human protein atlas: A spatial map of the human proteome.
#'Protein Sci. 2018 Jan;27(1):233-244. doi: 10.1002/pro.3307.
#'Epub 2017 Oct 10. PMID: 28940711; PMCID: PMC5734309.\url{ https://doi.org/10.1002/pro.3307}
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
"tissue_list_protein"



#' Default list of RNA tissues used for correlation analysis
#'
#'A character vector of the default tissues that function correlation_per_gene
#'uses to analyze RNA vs protein expression when the input is a list of genes.
#'These tissue names follow the format of gtexr and can only be used to query
#'RNA expression.
#'
#'@format A character vector of length 10 (one element per tissue)
#'
#'@source Internal dataset created for ExprCompareR package.
#'Data derived from GTEx Portal (\url{https://gtexportal.org/}) via gtexr
#'
#'@examples
#'data(tissue_list_RNA)
#'length(tissue_list_RNA)
#'head(tissue_list_RNA)
#'\dontrun{
#'tissue_list_RNA
#'}
#'
#'@references
#'
#'Lonsdale, J., Thomas, J., Salvatore, M. et al.
#' The Genotype-Tissue Expression (GTEx) project. Nat Genet 45, 580-585 (2013).
#'\url{https://doi.org/10.1038/ng.2653}
#'
#'Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). "gtexr: A
#'convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API."
#'Journal of Open Source So ware, 10(109), 8249. ISSN
#'2475-9066, \url{doi:10.21105/joss.08249}, gigs v0.2.1.
#'
"tissue_list_RNA"



#' Tissue mapping: HPAanalyze to gtexr
#'
#' A tibble that maps tissue names from the Human Protein Atlas (HPAanalyze)
#' to the corresponding tissue names used in the GTEx portal via the gtexr package.
#' This is useful for converting HPA tissue names to GTEx-compatible names
#' when querying RNA expression data.
#'
#'@format A tibble with 33 rows and 2 variables:
#' \describe{
#'   \item{protein_tissue}{Tissue name used in HPAanalyze (character).}
#'   \item{RNA_tissue}{Corresponding tissue name used in GTEx (character).}
#' }
#'
#'@source  Internal dataset created for ExprCompareR package.
#'Data from Human Protein Atlas (\url{https://www.proteinatlas.org/}) via HPAanalyze
#'  and GTEx Portal (\url{https://gtexportal.org/}) via gtexr
#'
#'@examples
#' data(tissue_map)
#' length(tissue_map)
#' head(tissue_map)
#' tissue_map$protein_tissue
#' tissue_map$RNA_tissue
#' \dontrun{
#' tissue_map
#' }
#'
#'@references
#'
#'Lonsdale, J., Thomas, J., Salvatore, M. et al.
#'The Genotype-Tissue Expression (GTEx) project. Nat Genet 45, 580-585 (2013).
#'\url{https://doi.org/10.1038/ng.2653}
#'
#'Thul PJ, Lindskog C. The human protein atlas: A spatial map of the human proteome.
#'Protein Sci. 2018 Jan;27(1):233-244. doi: 10.1002/pro.3307.
#'Epub 2017 Oct 10. PMID: 28940711; PMCID: PMC5734309. \url{https://doi.org/10.1002/pro.3307}
#'
#'Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).
#'"HPAanalyze: an R package that facilitates the retrieval and analysis of the
#'Human Protein Atlas data." MC Bioinformatics 20, 463 (2019).
#'\url{https://doi.org/10.1186/s12859-019-3059-z}
#'
#'Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). "gtexr: A
#'convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API."
#'Journal of Open Source So ware, 10(109), 8249. ISSN
#'2475-9066, \url{doi:10.21105/joss.08249}, gigs v0.2.1.
#'
#'
"tissue_map"

