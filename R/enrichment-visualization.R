#' @import data.table
#' @import ggplot2
#' @importFrom corrplot corrplot
#' @importFrom ComplexHeatmap Heatmap
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MA plot
###
.maplot <- function(dt, pvalue, alpha,
                    xlab, ylab, title, subtitle, caption) {

    p <- ggplot(dt,aes(x = A, y = M, colour = get(pvalue) < alpha)) +
        geom_point() +
        labs(x = xlab, y = ylab, colour = "significant",
             title = title, subtitle = subtitle, caption = caption)

    if ("subset" %in% names(dt))
        p <- p + facet_wrap( ~ subset)

    return(p)
}

#' @export
setGeneric(
    name = "enrichmaplot",
    def = function(x, variable, pvalue = "pvals", alpha = 0.05,
                   xlab = "A", ylab = "M",
                   title = character(), subtitle = character(), caption = character()) {
        standardGeneric("enrichmaplot")
    }
)
#' @export
setMethod(
    f = "enrichmaplot",
    signature = c(x = "list"),
    definition = function(x, pvalue, alpha,
                          xlab, ylab, title, subtitle, caption) {

        dt <-data.table::rbindlist(x, fill = TRUE, idcol = "subset")

        .maplot(dt, pvalue, alpha,
                 xlab, ylab, title, subtitle, caption)
    }
)
#' @export
setMethod(
    f = "enrichmaplot",
    signature = c(x = "data.frame"),
    definition = function(x, pvalue, alpha,
                          xlab, ylab, title, subtitle, caption) {

        if (!("data.table" %in% class(x)))
            x <- as.data.table(x)

        .maplot(x, pvalue, alpha,
                xlab, ylab, title, subtitle, caption)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Barplot
###

.barplot <- function(dt, variable, pvalue, alpha,
                     xlab, ylab, title, subtitle, caption,
                     xangle, xsize, vjust,
                     yangle, ysize, hjust) {

    if (length(alpha) == 1)
        dt <- dt[get(pvalue) < alpha]
    else if (length(alpha) > 1)
        stop("A single value of alpha should be specified!")

    #if (length(threshold) == 1)
    #    dt <- dt[get(variable) > threshold,]

    #if (length(percentiles) == 1)
    #    dt <- dt[(.N-percentiles*.N):.N]

    p <- ggplot(dt, aes(x = reorder(category, -get(variable)), y = get(variable), fill = get(pvalue))) +
        geom_bar(stat = "identity", position=position_dodge()) +
        labs(x = xlab, y = ylab, fill = pvalue,
             title = title, subtitle = subtitle, caption = caption) +
        theme(axis.text.x  = element_text(angle=xangle, size=xsize, vjust=vjust),
              axis.text.y  = element_text(angle=yangle, size=ysize, hjust=hjust))

    if ("subset" %in% names(dt))
        p <- p + facet_grid(subset ~ .)

    return(p)
}

#' @export
setGeneric(
    name = "enrichbarplot",
    def = function(x, variable, pvalue = "pvals", alpha = numeric(),
                   xlab = character(), ylab = character(),
                   title = character(), subtitle = character(), caption = character(),
                   xangle = 90, xsize = 10, vjust = 1,
                   yangle = 0, ysize = 10, hjust = 1) {
        standardGeneric("enrichbarplot")
    }
)
#' @export
setMethod(
    f = "enrichbarplot",
    signature = c(x = "list"),
    definition = function(x, variable, pvalue, alpha,
                          xlab, ylab, title, subtitle, caption,
                          xangle, xsize, vjust,
                          yangle, ysize, hjust) {

        out <- sapply(x, function(y){
            y[, c("category", variable, pvalue), with = FALSE]
        }, simplify = FALSE, USE.NAMES = TRUE)

        dt <-data.table::rbindlist(out, fill = TRUE, idcol = "subset")
        keycols <- c("subset", variable)
        setkeyv(dt, keycols)

        .barplot(dt, variable, pvalue, alpha,
                 xlab, ylab, title, subtitle, caption,
                 xangle, xsize, vjust,
                 yangle, ysize, hjust)
    }
)
#' @export
setMethod(
    f = "enrichbarplot",
    signature = c(x = "data.frame"),
    definition = function(x, variable, pvalue, alpha,
                          xlab, ylab, title, subtitle, caption,
                          xangle, xsize, vjust,
                          yangle, ysize, hjust) {

        if (!("data.table" %in% class(x)))
            x <- as.data.table(x)

        dt <- x[, c("category", variable, pvalue), with = FALSE]
        setkeyv(dt, variable)

        .barplot(dt, variable, pvalue, alpha,
                 xlab, ylab, title, subtitle, caption,
                 xangle, xsize, vjust,
                 yangle, ysize, hjust)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sample correlation plot
###


.makemat <- function(x, variable, replace.na) {
    out <- lapply(1:length(x), function(y){
        DT <- x[[y]][, c("category", variable), with = FALSE]
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

.pvalmat <- function(x, pvalue) {
    if (!is.null(pvalue)) {
        p.mat <- lapply(1:length(x), function(y){
            DT <- x[[y]][, c("category", pvalue), with = FALSE]
            setnames(DT, pvalue, names(x)[y])
        })
        p.mat <- Reduce(function(...) merge(..., all = TRUE), p.mat)
    } else {
        p.mat = NULL
    }
    return(p.mat)
}

# @param pvalue A character vector, either "pvals" or "padj".
#' @export
setGeneric(
    name = "correlationplot",
    def = function(x, variable, replace.na = FALSE, method = "color", add = FALSE, col = NULL, bg = "white", title = "",
                   outline = FALSE, mar = c(0, 0, 0, 0), addgrid.col = NULL,
                   tl.pos = "lt", tl.cex = 0.6, tl.col = "black", tl.offset = 0.4, tl.srt = 90,
                   cl.pos = "b", cl.lim = NULL, cl.length = 2, cl.cex = 0.6,
                   cl.ratio = 0.15, cl.align.text = "c", cl.offset = 0.5,
                   addshade = "all", shade.lwd = 1, shade.col = shade.col,
                   pvalue = NULL, sig.level = sig.level, insig = "n",
                   pch = 4, pch.col = "black", pch.cex = 3,
                   na.label = "square", na.label.col = "white"){
        standardGeneric("correlationplot")
    }
)
#' @export
setMethod(
    f = "correlationplot",
    signature = "list",
    definition = function(x, variable, replace.na,
                          method, add, col, bg, title,
                          outline, mar, addgrid.col,
                          tl.pos, tl.cex, tl.col, tl.offset, tl.srt,
                          cl.pos, cl.lim, cl.length, cl.cex,
                          cl.ratio, cl.align.text, cl.offset,
                          addshade, shade.lwd, shade.col,
                          pvalue, sig.level, insig,
                          pch, pch.col, pch.cex,
                          na.label, na.label.col) {
        # if nested list, unlist elements which are lists
        nl <- lapply(x, class) == "list"
        if (any(nl)) {
            x <- unlist(x[nl], recursive = F)
            x <- c(x, x[!nl])
        }

        #TO-DO check that elements of the list are data.tables w. enrichment results

        dm <- .makemat(x, variable)
        p.mat <- .pvalmat (x, pvalue)


        corrplot(dm, replace.na, method, type = "full", add = add, col = col, bg = bg, title = "",
                 is.corr = FALSE, diag = FALSE, outline = outline, mar = mar,
                 addgrid.col = addgrid.col, addCoef.col = NULL, order = "original",
                 tl.pos = tl.pos, tl.cex = tl.cex, tl.col = tl.col, tl.offset = tl.offset, tl.srt = tl.srt,
                 cl.pos = cl.pos,cl.lim = cl.lim, cl.length = cl.length, cl.cex = cl.cex,
                 cl.ratio = cl.ratio, cl.align.text = cl.align.text, cl.offset = cl.offset,
                 addshade = "all", shade.lwd = 1, shade.col = shade.col,
                 p.mat = p.mat, sig.level = sig.level, insig = insig,
                 pch = pch, pch.col = pch.col, pch.cex = pch.cex,
                 plotCI = "n", na.label = na.label, na.label.col = na.label.col)

    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Heatmap
###

#' @export
setGeneric(
    name = "enrichheatmap",
    def = function(x, variable, replace.na = TRUE, ...){
        standardGeneric("enrichheatmap")
    }
)
#' @export
setMethod(
    f = "enrichheatmap",
    signature = "list",
    definition = function(x, variable, ...) {

        # if nested list, unlist elements which are lists
        nl <- lapply(x, class) == "list"
        if (any(nl)) {
            x <- unlist(x[nl], recursive = F)
            x <- c(x, x[!nl])
        }
        #TO-DO check that elements of the list are data.tables w. enrichment results

        dm <- .makemat(x, variable, replace.na) # a matrix for plotting heatmap
        # p.mat <- .pvalmat (x, pvalue)
        Heatmap(matrix = dm, ...)
    }
)
