#' @include enrich.data.frame-class.R
#' @import data.table
#' @import ggplot2
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MA plot
###

.maplot <- function(dt, pvalue, alpha) {

    p <- ggplot(dt,aes(x = A, y = M, colour = get(pvalue) < alpha)) +
        geom_point() +
        labs(x = "A", y = "M", colour = "significant")

    if ("subset" %in% names(dt))
        p <- p + facet_wrap( ~ subset)

    return(p)
}

#' @export
setGeneric(
    name = "enrichMAplot",
    def = function(x, variable, pvalue = "pvals", alpha = 0.05) {
        standardGeneric("enrichMAplot")
    }
)
#' MA plot of enriched annotations.
#'
#' Make an MA-like plot of enriched annotations, similar to the commonly used
#' plots in differential expression analysis.
#'
#' @param x \code{enrich.data.frame} object, or a list of those.
#' @param pvalue Character, one of \code{c("pvals", "padj")}.
#' @param alpha Numeric, significance level to be used for plotting.
#'
#' @return A \code{ggplot} object.
#'
#' @name enrichMAplot
#' @export
setMethod(
    f = "enrichMAplot",
    signature = c(x = "list"),
    definition = function(x, pvalue, alpha) {

        if (!all(sapply(x, class) == "list"))
            stop("x should be a list of enrich.data.frame objects!")

        dt <- data.table::rbindlist(x, fill = TRUE, idcol = "subset")

        .maplot(dt, pvalue, alpha)
    }
)

#' @rdname enrichMAplot
#' @export
setMethod(
    f = "enrichMAplot",
    signature = c(x = "enrich.data.frame"),
    definition = function(x, pvalue, alpha) {

        .maplot(x, pvalue, alpha)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Barplot
###

.barplot <- function(dt, variable, pvalue, alpha) {

    if (length(alpha) == 1)
        dt <- dt[get(pvalue) < alpha]
    else if (length(alpha) > 1)
        stop("A single value of alpha should be specified!")

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

#' @export
setGeneric(
    name = "enrichBarplot",
    def = function(x, variable, pvalue = "pvals", alpha = numeric()) {
        standardGeneric("enrichBarplot")
    }
)
#' Barplot of enriched and depleted annotations.
#'
#' Make a barplot of enriched annotations. Bars' heights represent values of
#' the chosen enrichment statistic (\code{c("enrich","M","A")}), and the colours
#' represent the p values (\code{c("pvals", "padj")}).
#'
#' @inheritParams enrichMAplot
#' @param variable Character, indicating the statistic values to be used for
#'    plotting, must be one of \code{c("enrich","M","A")}.
#'
#' @return A \code{ggplot} object.
#'
#' @name enrichBarplot
#' @export
setMethod(
    f = "enrichBarplot",
    signature = c(x = "list"),
    definition = function(x, variable, pvalue, alpha) {

        if (!all(sapply(x, class) == "list"))
            stop("x should be a list of enrich.data.frame objects!")

        out <- sapply(x, function(y){
            y[, c("category", variable, pvalue), with = FALSE]
        }, simplify = FALSE, USE.NAMES = TRUE)

        dt <- data.table::rbindlist(out, fill = TRUE, idcol = "subset")
        keycols <- c("subset", variable)
        setkeyv(dt, keycols)

        .barplot(dt, variable, pvalue, alpha)
    }
)
#' @rdname enrichBarplot
#' @export
setMethod(
    f = "enrichBarplot",
    signature = c(x = "enrich.data.frame"),
    definition = function(x, variable, pvalue, alpha) {

        dt <- x[, c("category", variable, pvalue)]
        setkeyv(dt, variable)

        .barplot(dt, variable, pvalue, alpha)
    }
)

