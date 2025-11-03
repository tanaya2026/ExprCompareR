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
#'@references
#'
#'Lonsdale, J., Thomas, J., Salvatore, M. et al.
#'The Genotype-Tissue Expression (GTEx) project. Nat Genet 45, 580–585 (2013).
#'https://doi.org/10.1038/ng.2653
#'
#'Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). “gtexr: A
#'convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API.”
#'Journal of Open Source So ware, 10(109), 8249. ISSN
#'2475-9066, doi:10.21105/joss.08249, gigs v0.2.1.
#'
#'@examples
#'data(tissue_list_RNA)
#'length(tissue_list_RNA)
#'head(tissue_list_RNA)
#'\dontrun{
#'tissue_list_RNA
#'}
"tissue_list_RNA"
