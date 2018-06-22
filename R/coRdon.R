#' coRdon: codon usage analysis in R
#'
#' R package for analysis of codone usage in unannotated or
#' KEGG/COG annotated DNA sequences. Calculates various
#' measures of CU bias and CU-based predictors of gene
#' expression, and performs gene set enrichment analysis for
#' annotated sequences. Implements several methods for
#' visualization of CU and enrichment analysis results.
#'
#' @docType package
#' @name coRdon
#' @importFrom methods setClass setGeneric setMethod setValidity validObject
#' @importFrom methods show is new
#' @importFrom stats binom.test p.adjust quantile reorder
#' @importFrom utils capture.output head str tail unzip sessionInfo
NULL

# quiets NOTEs in R CMD check:
# column and row names for objects defined in the package
#' @importFrom utils globalVariables
globalVariables(c("category", "x", "y", "genes", "AA", "Start", "annot",
                "A", "M", "CATEGORY", "ANN", "codon", "..row"))
