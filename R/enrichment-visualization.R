#' @import data.table
#' @importFrom corrplot corrplot
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Barplot
###

.barplot <- function(dt, variable, pvals, threshold,
                     xlab, ylab, title, subtitle, caption,
                     xangle, xsize, vjust,
                     yangle, ysize, hjust) {

    if (length(threshold) == 1)
        dt <- dt[get(variable) > threshold,]

    #if (length(percentiles) == 1)
    #    dt <- dt[(.N-percentiles*.N):.N]

    ggplot(dt, aes(x = reorder(category, -get(variable)), y = get(variable), fill= pvals)) +
        geom_bar(stat = "identity", position=position_dodge()) +
        facet_grid(subset ~ .) +
        labs(x = xlab, y = ylab,
             title = title, subtitle = subtitle, caption = caption) +
        theme(axis.text.x  = element_text(angle=xangle, size=xsize, vjust=vjust),
              axis.text.y  = element_text(angle=yangle, size=ysize, hjust=hjust))
}

#' @export
setGeneric(
    name = "enrichbarplot",
    def = function(x, variable, pvals = NULL, threshold = numeric(),
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
    definition = function(x, variable, pvals, threshold,
                          xlab, ylab, title, subtitle, caption,
                          xangle, xsize, vjust,
                          yangle, ysize, hjust) {

        out <- sapply(x, function(y){
            y[, c("category", variable, pvals), with = FALSE]
        }, simplify = FALSE, USE.NAMES = TRUE)

        dt <-data.table::rbindlist(out, fill = TRUE, idcol = "subset")
        keycols <- c("subset", variable)
        setkeyv(dt, keycols)

        .barplot(dt, variable, pvals, threshold,
                 xlab, ylab, title, subtitle, caption,
                 xangle, xsize, vjust,
                 yangle, ysize, hjust)
    }
)
#' @export
setMethod(
    f = "enrichbarplot",
    signature = c(x = "data.frame"),
    definition = function(x, variable, pvals, threshold,
                          xlab, ylab, title, subtitle, caption,
                          xangle, xsize, vjust,
                          yangle, ysize, hjust) {

        if (!("data.table" %in% class(x)))
            x <- as.data.table(x)

        dt <- x[, subset := NA][, c("subset", "category", variable, pvals), with = FALSE]
        keycols <- c("subset", variable)
        setkeyv(dt, keycols)

        .barplot(dt, variable, pvals, threshold,
                 xlab, ylab, title, subtitle, caption,
                 xangle, xsize, vjust,
                 yangle, ysize, hjust)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sample correlation plot
###

# @param pvals A character vector, either "pvals" or "padj".
#' @export
setGeneric(
    name = "correlationplot",
    def = function(x, variable, method = "color", add = FALSE, col = NULL, bg = "white", title = "",
                   outline = FALSE, mar = c(0, 0, 0, 0), addgrid.col = NULL,
                   tl.pos = "lt", tl.cex = 0.6, tl.col = "black", tl.offset = 0.4, tl.srt = 90,
                   cl.pos = "b", cl.lim = NULL, cl.length = 2, cl.cex = 0.6,
                   cl.ratio = 0.15, cl.align.text = "c", cl.offset = 0.5,
                   addshade = "all", shade.lwd = 1, shade.col = shade.col,
                   pvals = NULL, sig.level = sig.level, insig = "n",
                   pch = 4, pch.col = "black", pch.cex = 3,
                   na.label = "square", na.label.col = "white"){
        standardGeneric("correlationplot")
    }
)
#' @export
setMethod(
    f = "correlationplot",
    signature = "list",
    definition = function(x, variable,
                          method, add, col, bg, title,
                          outline, mar, addgrid.col,
                          tl.pos, tl.cex, tl.col, tl.offset, tl.srt,
                          cl.pos, cl.lim, cl.length, cl.cex,
                          cl.ratio, cl.align.text, cl.offset,
                          addshade, shade.lwd, shade.col,
                          pvals, sig.level, insig,
                          pch, pch.col, pch.cex,
                          na.label, na.label.col) {
        # if nested list, unlist elements which are lists
        nl <- lapply(x, class) == "list"
        if (any(nl)) {
            x <- unlist(x[nl], recursive = F)
            x <- c(x, x[!nl])
        }

        #TO-DO check that elements of the list are data.tables w. enrichment results

        out <- lapply(1:length(x), function(y){
            DT <- x[[y]][, c("category", variable), with = FALSE]
            setnames(DT, variable, names(x)[y])
        })
        dt <- Reduce(function(...) merge(..., all = TRUE), out)
        # dt[is.na(dt)] <- 0
        dm <- data.matrix(dt[,-1])
        rownames(dm) <- unname(unlist(dt[,1]))

        if (!is.null(pvals)) {
            p.mat <- lapply(1:length(x), function(y){
                DT <- x[[y]][, c("category", pvals), with = FALSE]
                setnames(DT, pvals, names(x)[y])
            })
            p.mat <- Reduce(function(...) merge(..., all = TRUE), p.mat)
        } else {
            p.mat = NULL
        }

        corrplot(dm, method, type = "full", add = add, col = col, bg = bg, title = "",
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
