#' @include enrichment.R
NULL

#' An S4 class \code{crossTab}
#'
#' A contingency table that displays...
#'
#' @slot sequences A character vector of sequences annotations (KO, COG).
#' @slot variable A coresponding numeric vector of CU values.
#' @slot table A contingecy table.
#'
setClass(
    "crossTab",
    slots = c(
        sequences = "character",
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
    def = function(sequences, variable, ...){
        standardGeneric("crossTab")
    }
)

#' @describeIn crossTab Create new objects of class \code{crossTab}.
#'
#' Creates a contingency table for the list of annotated sequences and
#' the corresponding list of codon usage (CU) values.
#'
#' @param sequences A character vector of sequences annotations (KO, COG).
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
    signature = c(sequences = "character", variable = "numeric"),
    definition = function(sequences, variable, threshold = 1L, percentiles = NULL) {
        new("crossTab",
            sequences = sequences,
            variable = variable,
            table = make.contable(sequences, variable, threshold, percentiles))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### crossTab accesor methods
###

#' Display the object of \code{crossTab} class.
#'
#' @inheritParams methods::show
#'
#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = "crossTab",
    definition = function(object){
        l <- length(object@sequences)
        s <- capture.output(str(object@sequences))
        v <- capture.output(str(object@variable))
        tab <- object@table
        cat("crossTab instance for", l, "sequences.\n")
        cat("Sequence annotations:\n", s, "\n")
        cat("Corresponding variable values:\n", v, "\n")
        cat("Contingecy table:\n")
        print(tab)
    }
)

#' Length of \code{crossTab} object
#'
#' Get the length of \code{crossTab} object, i.e. the number of sequences
#' for which the contingency table is contained in the object.
#'
#' @param x A \code{crossTab} object.
#'
#' @rdname length
#' @export
setMethod(
    f = "length",
    signature = "crossTab",
    definition = function(x){
        length(x@sequences)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "getSeqAnnot",
    def = function(x){
        standardGeneric("getSeqAnnot")
    }
)

#' @describeIn crossTab Get sequence annotations from \code{crossTab} object.
#'
#' @param x A \code{crossTab} object.
#'
#' @export
setMethod(
    f = "getSeqAnnot",
    signature = "crossTab",
    definition = function(x){
        return(x@sequences)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "getVariable",
    def = function(x){
        standardGeneric("getVariable")
    }
)

#' @describeIn crossTab Get values of the variable used to create contingency table
#'    contained in \code{crossTab} object.
#'
#' @param x A \code{crossTab} object.
#'
#' @export
setMethod(
    f = "getVariable",
    signature = "crossTab",
    definition = function(x){
        return(x@variable)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "contable",
    def = function(x){
        standardGeneric("contable")
    }
)

#' @describeIn crossTab Get contingency table from \code{crossTab} object.
#'
#' @param x A \code{crossTab} object.
#'
#' @export
setMethod(
    f = "contable",
    signature = "crossTab",
    definition = function(x){
        return(x@table)
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
#' Reduce the input contingency table by associating sequences with KEGG Pathway
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
            sequences = x@sequences,
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
