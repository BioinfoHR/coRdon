#' @include enrich.data.frame-class.R
#' @import data.table
NULL

.enrichment <- function(contable, pvalueCutoff, pAdjustMethod, padjCutoff) {

    rows <- names(contable)
    top_rows <- rows[grep("top|gt", rows)]

    all <- contable$all
    all.sum <- sum(all)

    by.top <- sapply(top_rows, function(row) {

        top <- unname(unlist(contable[,..row]))
        top.sum <- sum(top)

        sc <- top.sum / all.sum

        scaled_top <- top + 1
        scaled_all <- all * sc + 1

        ct <- contable[, c("category", "all", row), with = FALSE]

        pvals <-
            apply(ct[,c("all", row), with = FALSE], 1, function(x) {
                b = binom.test(x[2], top.sum, x[1]/all.sum)
                b$p.value
            })
        padj <- p.adjust(pvals, method = pAdjustMethod)

        ct[, ':='(enrich = (scaled_top - scaled_all) / scaled_all * 100,
                  M = log2(scaled_top) - log2(scaled_all),
                  A = (log2(scaled_all) + log2(scaled_top)) / 2,
                  pvals = pvals,
                  padj = padj)]

        if(length(pvalueCutoff) != 0) ct <- ct[pvals <= pvalueCutoff,]
        if(length(padjCutoff) != 0) ct <- ct[padj <= padjCutoff,]

        enrich.data.frame(ct)

    }, simplify = FALSE, USE.NAMES = TRUE)
}

#' @export
setGeneric(
    name = "enrichment",
    def = function(x, pvalueCutoff = numeric(), pAdjustMethod = "BH", padjCutoff = numeric()){
        standardGeneric("enrichment")
    }
)

#' Enrichment analysis for codon usage (CU) data.
#'
#' Performs enrichment analysis, given a contongency table of codon counts.
#' p values are calculated by binomial test, adjustment for multiple testing
#' can be performed by any of the `p.adjust.methods`.
#'
#' @param x `crossTab` object
#' @param pvaluCutoff numeric, discard categories with p value below this
#'    threshold. By default, no threshold is set (`numeric()`).
#' @param pAdjustMethod character, one of the `p.adjust.methods`.
#' @param padjCutoff numeric, discard categories with adjusted p value below
#'    this threshold. By default, no threshold is set (`numeric()`).
#'
#' @return An `enrich.data.frame` object, or a list of those.
#'
setMethod(
    f = "enrichment",
    signature = "crossTab",
    definition = function(x, pvalueCutoff, pAdjustMethod, padjCutoff){
        enl <- .enrichment(x@table, pvalueCutoff, pAdjustMethod, padjCutoff)
        if (length(enl) == 1) enl[[1]]
        else enl
    }
)

.makemat <- function(x, variable, replace.na) {
    out <- lapply(1:length(x), function(y){
        DT <- x[[y]][, c("category", variable)]
        setnames(DT, variable, names(x)[y])
    })
    dt <- Reduce(function(...) merge(..., all = TRUE), out)
    if (replace.na) {
        if (is.logical(replace.na)) replace.na = 0
        for (j in seq_len(ncol(dt)))
            set(dt, which(is.na(dt[[j]])), j, replace.na)
    }
    dm <- data.matrix(dt[,-1])
    rownames(dm) <- unname(unlist(dt[,1]))
    return(dm)
}

#' Extract enrichment values from multiple samples. Note that the samples should
#' contain annotations of the same type (i.e. the same onthology).
#'
#' @param x list of `enrich.data.frame` objects
#' @param variable character string indicating the statistic values to extract
#'    from `enrich.data.frame` objects in x, must be one of `c("enrich","M","A")`.
#' @param replace.na logical, whether to replace NA values in the output.
#'    If `TRUE` (default), NAs will be replaced by 0. Alternatively, if numueric,
#'    NAs will be replaced by that given value.
#' @return `matrix` with sequences' annotations as rows, and variable values
#'   for different samples as columns.
#'
#' @export
setGeneric(
    name = "enrich.matrix",
    def = function(x, variable, replace.na = TRUE, ...){
        standardGeneric("enrich.matrix")
    }
)

#' @export
setMethod(
    f = "enrich.matrix",
    signature = c(x = "list"),
    definition = function(x, variable, replace.na){

        # if nested list, unlist elements which are lists
        nl <- lapply(x, class) == "list"
        if (any(nl)) {
            x <- unlist(x[nl], recursive = F)
            x <- c(x, x[!nl])
        }

        # enrich.data.frame class
        if (!(all(sapply(x, class) == "enrich.data.frame")))
            stop("x should be a (nested) list of enrich.data.frame objects!")

        .makemat(x, variable, replace.na)
    }
)
