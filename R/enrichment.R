#' @include codonTable-class.R
#' @include genCode-class.R
#' @import data.table
NULL

#' Make contingency table
#'
#' Creates a contingency table for the list of annotated sequences and
#' the corresponding list of codon usage (CU) values.
#'
#' @param geneList A character vector of sequences annotations (KO, COG).
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
#'
make.contable <- function(geneList, variable,
                          threshold = 1L, percentiles = NULL) {

  # geneList <- as.factor(slot(cTobject, category))
  geneList <- as.factor(geneList)
  all <- as.vector(table(geneList))
  result <- data.table(category = levels(geneList),
                       all = all)

  if(!is.null(percentiles)) {
    top.perc <- lapply(percentiles, function(x) {
      as.vector(table(geneList[variable >= quantile(variable, x)]))
    })
    names <- make.names(paste("top", percentiles, sep="_"))
    result[, (names) := top.perc]
  } else {
    top.perc <- NULL
  }

  if(!is.null(threshold)) {
    top.thresh <- lapply(threshold, function(x) {
      as.vector(table(geneList[variable >= x]))
    })
    names <- make.names(paste("gt", threshold, sep="_"))
    result[, (names) := top.thresh]
  } else {
    top.thresh <- NULL
  }

  return(result)
}

#' Reduce contingency table
#'
#' Given contingency table with gene identifiers, reduce the table by associating
#' gene identifiers with either KEGG Pathway or KEGG Module identifiers.
#'
#' @param contable A contingency table.
#' @param target A character vector indicating which onthology to use, either
#'    \code{"pathway"} or \code{"module"}.
#'
#' @return Returns a data.table with new category values in rows, and with the same
#'    columns as in the input contingency table.
#'
#' @export
#'
reduce.contable <- function(contable, target) {
    if (target == "pathway") {
        DT <- KO_PATHWAYS
    } else if (target == "module") {
        DT <- KO_MODULES
    }
    values <- unique(DT[,CATEGORY])
    tt <- sapply(values, function(x){
        KOs <- DT[CATEGORY == x, KO]
        if (any(KOs %in% contable[,category]))
            contable[category %in% KOs, lapply(.SD, sum), .SDcols = names(contable)[-1]]
        else NULL
    }, simplify = FALSE, USE.NAMES = TRUE)
    out <- Filter(Negate(is.null), tt)
    rbindlist(out, idcol = "category")
}

#' Enrichment analysis
#'
#' Enrichment analysis...
#'
#' @param contable A contingency table.
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
enrichment <- function(contable) {

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
        padj <- p.adjust(pvals, method = "BH")

        ct[, ':='(enrich = (scaled_top - scaled_all) / scaled_all * 100,
                  M = log2(scaled_top) - log2(scaled_all),
                  A = (log2(scaled_all) + log2(scaled_top)) / 2,
                  pvals = pvals,
                  padj = padj)]

    }, simplify = FALSE)

    #bound <- do.call(cbind, by.top)
    #out <- data.frame(all = all, bound)

}
