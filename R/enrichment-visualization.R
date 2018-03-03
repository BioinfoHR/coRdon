#' @import data.table
#' @importFrom corrplot corrplot
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### correlation plot
###

# @param pvals A character vector, either "pvals" or "padj".
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
