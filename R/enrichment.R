#' @include codonTable-class.R
#' @include genCode-class.R
#' @import data.table
NULL

make.contable <- function(genes, variable,
                      threshold = 1L, percentiles = NULL) {

  # genes <- as.factor(slot(cTobject, category))
  genes <- as.factor(genes)
  all <- as.vector(table(genes))
  result <- data.table(category = levels(genes),
                       all = all)

  if(!is.null(percentiles)) {
    top.perc <- lapply(percentiles, function(x) {
      as.vector(table(genes[variable >= quantile(variable, x)]))
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


# #' @importClassesFrom DOSE enrichResult
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
        ct

    }, simplify = FALSE)
}
