#' Converts protein expression data queried from HPAanalyze, to numeric values
#'
#' A function that converts the protein expression data of a given gene/tissue
#' from HPAanalyze, into numeric values. HPAnalyze returns protein expression as
#' strings, for e.g.("low", "medium", etc). In order to determine the pearson
#' correlation between the RNA expression and protein expression, we convert these
#' strings into a numeric value.
#'
#' @param name description
#' @return
#' @examples
#' @references
#' @import
#'
