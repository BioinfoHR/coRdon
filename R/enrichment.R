#' @import data.table
#' @importClassesFrom Biobase AnnotatedDataFrame
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom Biobase pData
NULL

.enrichment <- function(contable, pvalueCutoff, pAdjustMethod, padjCutoff) {

    rows <- names(contable)
    top_rows <- rows[grep("top|gt", rows)]

    all <- contable$all
    all.sum <- sum(all)

    by.top <- lapply(top_rows, function(row) {

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

        df <- data.frame(labelDescription = names(ct))
        AnnotatedDataFrame(data = ct, varMetadata = df)

    })
    names(by.top) <- top_rows
    by.top
}

#' Enrichment analysis for codon usage (CU) data.
#'
#' Performs enrichment analysis, given a contongency table of codon counts.
#' p values are calculated by binomial test, adjustment for multiple testing
#' can be performed by any of the \code{p.adjust.methods}.
#'
#' @param x A \code{crossTab} object
#' @param pvalueCutoff Numeric, discard categories with p value below this
#'    threshold. By default, no threshold is set (\code{numeric()}).
#' @param pAdjustMethod Character, one of the \code{p.adjust.methods}.
#' @param padjCutoff Numeric, discard categories with adjusted p value below
#'    this threshold. By default, no threshold is set (\code{numeric()}).
#'
#' @return An \code{AnnotatedDataFrame} object, or a list of those; data in
#' each object has category values in rows, and the following columns:
#'    \itemize{
#'      \item category, a character vector of annotation categories
#'      \item all, a numeric vector of integers, coresponding to sequence
#'          counts for each annotation category, in the background gene set
#'          (universe).
#'      \item a numeric vector(s) of integers, coresponding to sequence counts
#'          for each annotation category, in the set of genes for which
#'          enrichment is calculated, i.e. the predefined subset of (usually
#'          highly expressed) genes in the universe (named for the
#'          corresponding `crossTab` column).
#'      \item enrichment, calculated as the ratio: (scaled sample counts -
#'          scaled backg. counts) / scaled backg. counts * 100,
#'          where scaling means that sample counts are simply increased by 1,
#'          and background counts are multiplied by ratio of summed sample
#'          counts and summed backgroun counts, and also increased by 1
#'      \item M, log ratios of scaled counts
#'      \item A, mean average of scaled counts
#'      \item pvals, p values for exact binomial test
#'      \item padj, p values corrected by BH method.
#'    }
#'
#' @examples
#' require(Biobase)
#'
#' # create contingency table
#' s <- getKO(HD59)
#' v <- as.numeric(MELP(HD59, ribosomal = TRUE))
#' ct <- crossTab(s, v)
#'
#' # enrichment analysis
#' enr <- enrichment(ct)
#' enr # for help, see `?Biobase::AnnotatedDataFrame`
#' head(pData(enr))
#'
#' enr <- enrichment(ct, pAdjustMethod = "holm")
#' head(pData(enr))
#'
#' enr <- enrichment(ct, pvalueCutoff = 0.05)
#' head(pData(enr))
#'
#' enr <- enrichment(ct, padjCutoff = 0.05)
#' head(pData(enr))
#'
#' @rdname enrichment
#' @export
#' @export
setGeneric(
    name = "enrichment",
    def = function(x, pvalueCutoff = numeric(),
                   pAdjustMethod = "BH", padjCutoff = numeric()){
        standardGeneric("enrichment")
    }
)

#' @rdname enrichment
#' @export
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
    out <- lapply(seq_along(x), function(y){
        DT <- as.data.table(pData(x[[y]][, c("category", variable)]))
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

#' Extract chosen enrichment values to a matrix.
#'
#' Extract enrichment values from multiple samples, i.e.
#' \code{AnnotatedDataFrame} objects. Note that the samples should contain
#' annotations of the same type (i.e. the same ontology). The data in matrix
#' format can be easily used in different types of downstream analyses,
#' such as GAGE, and visualised, e.g. using a heatmap.
#'
#' @param x A named list of \code{AnnotatedDataFrame} objects.
#' @param variable Character, indicating the statistic values to extract from
#'    \code{enrich.data.frame} objects in x, must be one of
#'    \code{c("enrich","M","A")}.
#' @param replace.na logical, whether to replace NA values in the output.
#'    If `TRUE` (default), NAs will be replaced by 0. Alternatively,
#'    if numueric, NAs will be replaced by that given value.
#'
#' @return \code{matrix} with sequences' annotations as rows, and variable
#'   values for different samples as columns.
#'
#' @examples
#' require(Biobase)
#'
#' # create contingency table
#' s <- getKO(LD94)
#' v <- as.numeric(MELP(LD94, ribosomal = TRUE))
#' ct <- crossTab(s, v, percentiles = 0.2)
#'
#' # enrichment analysis
#' enr <- enrichment(ct)
#' enr # for help, see `?Biobase::AnnotatedDataFrame`
#' head(pData(enr$top_0.2), 10)
#' head(pData(enr$gt_1), 10)
#' enrm <- enrich.matrix(enr, "M")
#' head(enrm)
#'
#' @rdname enrich.matrix
#' @export
setGeneric(
    name = "enrich.matrix",
    def = function(x, variable, replace.na = TRUE){
        standardGeneric("enrich.matrix")
    }
)

#' @rdname enrich.matrix
#' @export
setMethod(
    f = "enrich.matrix",
    signature = c(x = "list"),
    definition = function(x, variable, replace.na){

        # if nested list, unlist elements which are lists
        nl <- lapply(x, class) == "list"
        if (any(nl)) {
            x <- unlist(x[nl], recursive = FALSE)
            x <- c(x, x[!nl])
        }

        # AnnotatedDataFrame class
        if (!(all(vapply(x, class,
                         character(length = 1)) == "AnnotatedDataFrame")))
            stop("x should be a (nested) list of AnnotatedDataFrame objects!")

        .makemat(x, variable, replace.na)
    }
)
