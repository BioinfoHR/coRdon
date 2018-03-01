#' @include enrichment.R
NULL

#' An S4 class \code{crossTab}
#'
#' A contingency table that displays...
#'
#' @slot genes A character vector of sequences annotations (KO, COG).
#' @slot variable A coresponding numeric vector of CU values.
#' @slot table A contingecy table.
#'
setClass(
    "crossTab",
    slots = c(
        genes = "character",
        variable = "numeric",
        table = "data.table"
    )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### crossTab constructor
###

#' @rdname crossTab-class
#' @export
setGeneric(
    name = "crossTab",
    def = function(genes, variable, ...){
        standardGeneric("crossTab")
    }
)

#' @describeIn crossTab Create new objects of class \code{crossTab}.
#'
#' Creates a contingency table for the list of annotated sequences and
#' the corresponding list of codon usage (CU) values.
#'
#' @param genes A character vector of sequences annotations (KO, COG).
#' @param variable A coresponding numeric vector of CU values.
#' @param threshold A threshold value (or a vector of values) of the variable.
#'    Sequences with value of the given variable greater than threshold are
#'    counted as a subset. Default is 1.
#' @param percentiles A single value or a vector of values between 0 and 1.
#'    Sequences with value of the given variable in the top percentiles are
#'    counted as a subset.
#'
#' @return Returns a data.table with category values in rows, and with separate
#'    columns for counts in background (all) and subsets, i.e. for diferrent
#'    thresholds/percentiles provided.
#'
#' @export
setMethod(
    f = "crossTab",
    signature = c(genes = "character", variable = "numeric"),
    definition = function(genes, variable, threshold = 1L, percentiles = NULL) {
        new("crossTab",
            genes = genes,
            variable = variable,
            table = make.contable(genes, variable, threshold, percentiles))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reduce crossTab
###

#' Reduce objects of \code{crossTab} class.
#'
#' @export
setGeneric(
    name = "reduceCrossTab",
    def = function(x, target){
        standardGeneric("reduceCrossTab")
    }
)
#' @describeIn reduceCrossTab Reduce \code{crossTab}.
#'
#' Reduce the input contingency table by associating genes with KEGG Pathway
#' or KEGG Module identifiers.
#'
#' @param x A \code{crossTab} object to be reduced.
#' @param target A character vector indicating which onthology to use, either
#'    \code{"pathway"} or \code{"module"}.
#'
#' @return Returns input \code{crossTab} object, with updated contingency table,
#'    displaying new category values in rows, and updated counts in columns.
#'
#' @export
#'
setMethod(
    f = "reduceCrossTab",
    signature = c("crossTab", "character"),
    definition = function(x, target){
        new("crossTab",
            genes = x@genes,
            variable = x@variable,
            table = reduce.contable(x@table, target))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Enrichment analysis
###

#' Enrichment analysis for object \code{crossTab} class.
#' @export
setGeneric(
    name = "enrichment",
    def = function(x, pvalueCutoff = numeric(), pAdjustMethod = "BH", padjCutoff = numeric()){
        standardGeneric("enrichment")
    }
)
#' @describeIn enrichment Enrichment analysis.
#'
#' @param x A \code{crossTab} object.
#' @param target A character vector indicating which onthology to use, either
#'    \code{"pathway"} or \code{"module"}.
#'
#' @return  Returns a list of data.tables, each containing categories in rows,
#'    and the following columns:
#'    \itemize{
#'      \item background counts
#'      \item counts for a given subset
#'      \item enrichment, calculated as the ratio:
#'          (scaled sample counts - scaled backg. counts) / scaled backg. counts * 100,
#'          where scaling means that sample counts are simply increased by 1,
#'          and background counts are multiplied by ratio of summed sample counts
#'          and summed backgroun counts, and also increased by 1.
#'      \item M, log ratios of scaled counts
#'      \item A, mean average of scaled counts
#'      \item pvals, p values for exact binomial test
#'      \item padj, p values corrected by BH method.
#'    }
#'
#' @export
#'
setMethod(
    f = "enrichment",
    signature = "crossTab",
    definition = function(x, pvalueCutoff, pAdjustMethod, padjCutoff){
        .enrichment(x@table, pvalueCutoff, pAdjustMethod, padjCutoff)
    }
)
