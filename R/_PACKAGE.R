#' ExprCompareR: Statistical comparison and visualization of transcript-protein relationships
#'
#'@description
#'`ExprCompareR` is an R package designed to streamline the exploration and quantification of the relationship
#'between RNA expression and protein expression across human tissues and conditions.
#'The package provides functions for correlation analysis, outlier detection, analyzing differential expression between diseased and healthy tissues, and plotting.
#'
#'@name ExprCompareR
#'@keywords internal
#'@author Tanaya Datar
#'
#' @seealso
#' \code{\link{compute_correlation}}, \code{\link{detect_outliers}}, \code{\link{correlation_genes_only}},
#' \code{\link{correlation_tissues_only}},\code{\link{correlation_genes_tissues}},
#' \code{\link{per_gene_plot}},\code{\link{per_tissue_plot}},\code{\link{compareCancerProtein}}
#'
#'@references
#'* BioRender.com. BioRender. Available at: \url{https://www.biorender.com} (accessed 2 November 2025).
#'* Cetinkaya-Rundel M,Cheng J, Grolemund G (2017).Customize your UI with HTML. \url{https://shiny.posit.co/r/articles/build/html-tags/}
#'* Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Aden-Buie G, Xie Y,Allen J, McPherson J, Dipert A, Borges B (2025).
#'shiny: Web Application Framework for R. R package version 1.11.1, \url{https://CRAN.R-project.org/package=shiny}
#'* Lonsdale, J., Thomas, J., Salvatore, M. et al. The Genotype-Tissue Expression (GTEx) project. Nat Genet 45, 580–585 (2013). \url{https://doi.org/10.1038/ng.2653}
#'* Müller K, Wickham H (2025). tibble: Simple Data Frames. R package version 3.3.0, \url{https://tibble.tidyverse.org/.}
#'* OpenAI.(2025). ChatGPT (GPT-5). Retrieved November 2, 2025, from \url{https://chatgpt.com/}
#'* R Core Team (2024). R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. \url{https://www.R-project.org/}
#'* Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2019).A multivariate Poisson-log normal mixture model for clustering transcriptome sequencing data.BMC Bioinformatics. 2019;20(1):394. \url{https://pubmed.ncbi.nlm.nih.gov/31311497/}
#'* Thul, P. J., & Lindskog, C. (2018). The human protein atlas: A spatial map of the human proteome. Protein science : a publication of the Protein Society, 27(1), 233–244. \url{https://doi.org/10.1002/pro.3307}
#'* Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019). “HPAanalyze: an R package that facilitates the retrieval and analysis of the Human Protein Atlas data.”MC Bioinformatics 20, 463 (2019). \url{https://doi.org/10.1186/s12859-019-3059-z}
#'* Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). “gtexr: A convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API.” Journal of Open Source Software, 10(109), 8249. ISSN 2475-9066, \url{doi:10.21105/joss.08249}, gigs v0.2.1.
#'* Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3 319-24277-4, \url{https://ggplot2.tidyverse.org.}
#'* Wickham H, François R, Henry L, Müller K, Vaughan D (2025). dplyr: A Grammar of Data Manipulation. R package version 1.1.4, \url{https://dplyr.tidyverse.org.}
#'* Wickham H, Henry L (2025). purrr: Functional Programming Tools. R package version 1.1.0, \url{https://purrr.tidyverse.org/}
#'* Wickham H, Vaughan D, Girlich M (2025). tidyr: Tidy Messy Data. R package version 1.3.1, \url{https://tidyr.tidyverse.org.}
#'
#'
#'
#'

NULL
