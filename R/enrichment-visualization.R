#' @importClassesFrom Biobase AnnotatedDataFrame
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom Biobase pData
#' @import data.table
#' @import ggplot2
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MA plot
###

.maplot <- function(dt, pvalue, siglev) {

    p <- ggplot(dt,aes(x = A, y = M, colour = get(pvalue) < siglev)) +
        geom_point() +
        labs(x = "A", y = "M", colour = "significant")

    if ("subset" %in% names(dt))
        p <- p + facet_wrap( ~ subset)

    return(p)
}

#' MA plot of enriched annotations.
#'
#' Make an MA-like plot of enriched annotations, similar to
#' the commonly used plots in differential expression analysis.
#'
#' @param x \code{AnnotatedDataFrame} object, or a list of those.
#' @param pvalue Character, one of \code{c("pvals", "padj")}.
#' @param siglev Numeric, significance level to be used for plotting.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' require(ggplot2)
#'
#' HD59_KO
#' enrichMAplot(HD59_KO)
#' enrichMAplot(HD59_KO, pvalue = "padj")
#' enrichMAplot(HD59_KO, siglev = 0.01)
#' enrichMAplot(HD59_KO, pvalue = "padj", siglev = 0.01)
#'
#' x <- list(disease = LD94_KO, healthy = HD59_KO)
#' enrichMAplot(x)
#'
#' @rdname enrichMAplot
#' @export
setGeneric(
    name = "enrichMAplot",
    def = function(x, pvalue = "pvals", siglev = 0.05) {
        standardGeneric("enrichMAplot")
    }
)

#' @rdname enrichMAplot
#' @export
setMethod(
    f = "enrichMAplot",
    signature = c(x = "list"),
    definition = function(x, pvalue, siglev) {

        if (!all(vapply(x, class,
                        character(length = 1)) == "AnnotatedDataFrame"))
            stop("x should be a list of AnnotatedDataFrame objects!")

        x <- lapply(x, function(x) as.data.table(pData(x)))
        dt <- data.table::rbindlist(x, fill = TRUE, idcol = "subset")

        .maplot(dt, pvalue, siglev)
    }
)

#' @rdname enrichMAplot
#' @export
setMethod(
    f = "enrichMAplot",
    signature = c(x = "AnnotatedDataFrame"),
    definition = function(x, pvalue, siglev) {

        .maplot(pData(x), pvalue, siglev)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Barplot
###

.barplot <- function(dt, variable, pvalue, siglev) {

    if (length(siglev) == 1)
        dt <- dt[get(pvalue) < siglev]
    else if (length(siglev) > 1)
        stop("A single value of siglev should be specified!")

    p <- ggplot(dt, aes(x = reorder(category, -get(variable)),
                        y = get(variable),
                        fill = get(pvalue))) +
        geom_bar(stat = "identity", position=position_dodge()) +
        labs(x = "category", y = variable, fill = pvalue) +
        theme(axis.text.x  = element_text(angle=90, size=10, vjust=1),
              axis.text.y  = element_text(angle=0, size=10, hjust=1))

    if ("subset" %in% names(dt))
        p <- p + facet_grid(subset ~ .)

    return(p)
}

#' Barplot of enriched and depleted annotations.
#'
#' Make a barplot of enriched annotations. Bars' heights represent values of
#' the chosen enrichment statistic (\code{c("enrich","M","A")}), and the
#' colours represent the p values (\code{c("pvals", "padj")}).
#'
#' @inheritParams enrichMAplot
#' @param variable Character, indicating the statistic values to be used for
#'    plotting, must be one of \code{c("enrich","M","A")}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' require(ggplot2)
#'
#' HD59_PATHWAYS
#' enrichBarplot(HD59_PATHWAYS, variable = "M",
#'               pvalue = "padj", siglev = 0.01) +
#'    labs(y = "pathway count\nlog ratios", x = "KEGG Pathway")
#'
#' x <- list(disease = LD94_PATHWAYS, healthy = HD59_PATHWAYS)
#' enrichBarplot(x, variable = "enrich", pvalue = "padj", siglev = 0.01) +
#'     labs(y = "relative enrichment", x = "KEGG Pathway")
#'
#' @rdname enrichBarplot
#' @export
setGeneric(
    name = "enrichBarplot",
    def = function(x, variable, pvalue = "pvals", siglev = numeric()) {
        standardGeneric("enrichBarplot")
    }
)

#' @rdname enrichBarplot
#' @export
setMethod(
    f = "enrichBarplot",
    signature = c(x = "list"),
    definition = function(x, variable, pvalue, siglev) {

        if (!all(vapply(x, class,
                        character(length = 1)) == "AnnotatedDataFrame"))
            stop("x should be a list of AnnotatedDataFrame objects!")

        out <- lapply(x, function(y){
            y <- as.data.table(pData(y))
            y[, c("category", variable, pvalue), with = FALSE]
        })

        dt <- data.table::rbindlist(out, fill = TRUE, idcol = "subset")
        keycols <- c("subset", variable)
        setkeyv(dt, keycols)

        .barplot(dt, variable, pvalue, siglev)
    }
)
#' @rdname enrichBarplot
#' @export
setMethod(
    f = "enrichBarplot",
    signature = c(x = "AnnotatedDataFrame"),
    definition = function(x, variable, pvalue, siglev) {

        dt <- as.data.table(pData(x))
        dt <- dt[, c("category", variable, pvalue), with = FALSE]
        setkeyv(dt, variable)

        .barplot(dt, variable, pvalue, siglev)
    }
)

