#' An S4 class \code{enrich.data.frame}
#'
#' CU based enrichment analysis results for a set of genes.
#' A data.table with category values in rows, and with separate columns
#' for the following statistics:
#'    \itemize{
#'      \item category, a character vector of annotation categories
#'      \item all, a numeric vector of integers, coresponding to sequence counts
#'          for each annotation category, in the background gene set (universe).
#'      \item a numeric vector of integers, coresponding to sequence counts
#'          for each annotation category, in the set of genes for which enrichment
#'          is calculated, i.e. the predefined subset of (usually highly expressed)
#'          genes in the universe (name for the corresponding `crossTab` column).
#'      \item enrichment, calculated as the ratio:
#'          (scaled sample counts - scaled backg. counts) / scaled backg. counts * 100,
#'          where scaling means that sample counts are simply increased by 1,
#'          and background counts are multiplied by ratio of summed sample counts
#'          and summed backgroun counts, and also increased by 1
#'      \item M, log ratios of scaled counts
#'      \item A, mean average of scaled counts
#'      \item pvals, p values for exact binomial test
#'      \item padj, p values corrected by BH method.
#'    }
#'
setClass(
    "enrich.data.frame",
    contains = "data.frame" # data.table issue #1881
)

setGeneric("enrich.data.frame",
           def = function(x) {
               standardGeneric("enrich.data.frame")
           })
setMethod(
    f = "enrich.data.frame",
    signature = "data.frame",
    definition = function(x) {
        new("enrich.data.frame", x)
    }
)
