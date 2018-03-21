#' @include codonTable-class.R
#' @include genCode-class.R
#' @import data.table
NULL

#' An S4 class \code{crossTab}
#'
#' Contingency table of sequences' annotations and the corresponding numeric values.
#'
#' @slot sequences Character vector of sequences annotations (KO, COG).
#' @slot variable Numeric vector of the coresponding CU values.
#' @slot table Contingecy table.
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

make.contable <- function(genes, variable,
                          threshold = 1L, percentiles = NULL) {

    # genes <- as.factor(slot(cTobject, category))
    genes <- as.factor(genes)
    all <- as.vector(table(genes))
    result <- data.table(category = levels(genes),
                         all = all)

    if(!is.null(percentiles)) {
        top.perc <- lapply(percentiles, function(x) {
            as.vector(table(genes[variable >= quantile(variable, 1-x)]))
        })
        names <- make.names(paste("top", percentiles, sep="_"))
        result[, (names) := top.perc]
    } else {
        top.perc <- NULL
    }

    if(!is.null(threshold)) {
        top.thresh <- lapply(threshold, function(x) {
            as.vector(table(genes[variable >= x]))
        })
        names <- make.names(paste("gt", threshold, sep="_"))
        result[, (names) := top.thresh]
    } else {
        top.thresh <- NULL
    }

    return(result)
}

#' @rdname crossTab-class
#' @export
setGeneric(
    name = "crossTab",
    def = function(sequences, variable, threshold = 1L, percentiles = NULL){
        standardGeneric("crossTab")
    }
)
#' @describeIn crossTab Create a contingency table for the set of annotated
#' sequences and the corresponding codon usage (CU) values.
#'
#' @param sequences Character vector of sequences' annotations (KO, COG).
#' @param variable Numeric vector of the coresponding CU values.
#' @param threshold A threshold value (or a vector of values) of the variable.
#'    Sequences with value of the given variable greater than threshold are
#'    taken as a subset. Default is 1. If no threshold should be set, specify
#'    \code{threshold = NULL}
#' @param percentiles A single value or a vector of values between 0 and 1.
#'    Sequences with value of the given variable in the top percentiles are
#'    taken as a subset. If no percentiles should be specified, the argument
#'    takes the value \code{NULL}.
#'
#' @return Returns a \code{crossTab} object with category values in rows, and
#'    with separate columns for counts in background (all) and subsets, i.e.
#'    for diferrent thresholds/percentiles provided.
#'
#' @export
setMethod(
    f = "crossTab",
    signature = c(sequences = "character", variable = "numeric"),
    definition = function(sequences, variable, threshold, percentiles) {
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
#' @param object A \code{crossTab} object.
#'
#' @name crossTab-show
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

#' Length of \code{crossTab} object, i.e. the number of
#' sequences for which the contingency table is contained in the object.
#'
#' @param x A \code{crossTab} object.
#'
#' @name crossTab-length
#' @export
setMethod(
    f = "length",
    signature = "crossTab",
    definition = function(x){
        length(x@sequences)
    }
)

#' @rdname crossTab-class
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

#' @rdname crossTab-class
#' @export
setGeneric(
    name = "getVariable",
    def = function(x){
        standardGeneric("getVariable")
    }
)

#' @describeIn crossTab Get values of the variable used to create contingency
#' table in \code{crossTab} object.
#'
#' @inheritParams getSeqAnnot
#'
#' @export
setMethod(
    f = "getVariable",
    signature = "crossTab",
    definition = function(x){
        return(x@variable)
    }
)

#' @rdname crossTab-class
#' @export
setGeneric(
    name = "contable",
    def = function(x){
        standardGeneric("contable")
    }
)

#' @describeIn crossTab Get contingency table from \code{crossTab} object.
#'
#' @inheritParams getSeqAnnot
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

reduce.contable <- function(contable, target) {
    if (target == "pathway") {
        DT <- KO_PATHWAYS
    } else if (target == "module") {
        DT <- KO_MODULES
    } else if (target == "cogfunction") {
        DT <- COGs
    }
    values <- unique(DT[,CATEGORY])
    tt <- sapply(values, function(x){
        anns <- DT[CATEGORY == x, ANN]
        if (any(anns %in% contable[,category]))
            contable[category %in% anns, lapply(.SD, sum), .SDcols = names(contable)[-1]]
        else NULL
    }, simplify = FALSE, USE.NAMES = TRUE)
    out <- Filter(Negate(is.null), tt)
    rbindlist(out, idcol = "category")
}

#' Reduce \code{crossTab}.
#'
#' Reduce the input contingency table by associating sequences with KEGG Pathway
#' or KEGG Module identifiers.
#'
#' @param x A \code{crossTab} object to be reduced.
#' @param target Character vector indicating which onthology to use, either
#'    \code{"pathway"} or \code{"module"}.
#'
#' @return Returns input \code{crossTab} object, with updated contingency table,
#'    displaying new category values in rows, and updated counts in columns.
#'
#' @rdname reduceCrossTab
#' @export
setGeneric(
    name = "reduceCrossTab",
    def = function(x, target){
        standardGeneric("reduceCrossTab")
    }
)

#' @rdname reduceCrossTab
#' @export
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

