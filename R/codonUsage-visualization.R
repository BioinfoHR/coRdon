#' @import data.table
#' @import ggplot2
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### B plot
###

.makedt <- function(dt, annotations, ribosomal, reference) {

    la <- length(annotations)
    if (la != 0){
        if (la != nrow(dt))
            stop(cat("Length of annotations vector", la,
                     "differs from the number of sequences,", nrow(dt)))
        else dt[, annotation := annotations]
    }

    dt[, genes := "other"]

    if (ribosomal){
        if (la == 0)
            stop("Ribosomal is TRUE, but no annotations were given!")
        else dt[annotation %in% RPKOs, genes := "ribosomal"]
    }

    if (length(reference) != 0){
        if (is.null(names(reference)))
            stop("Reference must be a named list!")
        if (all(unlist(lapply(reference, function(x) class(x) == "logical"))))
            invisible(
                lapply(names(reference), function(x){
                    dt[reference[[x]], genes := x]
                })
            )
        else if (all(unlist(lapply(reference, function(x) class(x) == "character"))))
            invisible(
                lapply(names(reference), function(x){
                    dt[annotation %in% reference[[x]], genes := x]
                })
            )
    }
    dt
}
.bplot <- function(dt){
    ggplot(dt, aes(x, y, colour = genes, alpha = genes)) +
        geom_jitter() +
        theme_light()
}

#' Karlin B plot
#'
#' Plot distances of each gene's CU frequency to specified gene (sub)sets (given by
#' \code{x} and \code{y}).
#'
#' @param x,y Character, both must be in \code{colnames(data)}, or numeric vectors
#'   of CU statistic values for two subsets of genes. If numeric, the vectors
#'   must be of the same length.
#' @param data A matrix with CU statistic values for subsets of genes in columns.
#' @param annotations A character vector giving KO annotations for sequences
#'   for which the CU values were calculated, must be of length \code{nrow(data)}.
#' @param ribosomal Logical, whether to indicate ribosomal genes in the plot.
#'   Default is \code{FALSE}, if set to \code{TRUE}, then \code{annotation} must
#'   be given.
#' @param reference A named list of logical vector(s) (each of length
#'   \code{nrow(data)}) of reference genes to be indicated on the plot, or a named
#'   list of character vector(s) (of any length) of the reference genes' anotations.
#'   If latter is the case, then \code{annotation} must be given.
#'
#' @return A \code{ggplot} object.
#'
#' @rdname bplot
#' @export
setGeneric(
    name = "bplot",
    def = function(x, y, data, annotations = character(),
                   ribosomal = FALSE, reference = list()){
        standardGeneric("bplot")
    }
)

#' @rdname bplot
#' @export
setMethod(
    f = "bplot",
    signature = c(x = "character", y = "character", data = "matrix"),
    definition = function(x, y, data, annotations, ribosomal, reference){

        dt <- as.data.table(data[, c(x,y)])

        if (length(xlab) == 0) xlab <- x
        if (length(ylab) == 0) ylab <- y

        setnames(dt, c(x, y), c("x","y"))
        dt <- .makedt(dt, annotations, ribosomal, reference)
        .bplot(dt)
    }
)
#' @rdname bplot
#' @export
setMethod(
    f = "bplot",
    signature = c(x = "numeric", y = "numeric", data = "missing"),
    definition = function(x, y, data, annotations, ribosomal, reference){

        dt <- data.table(cbind(x, y))

        if (length(xlab) == 0) xlab <- "x"
        if (length(ylab) == 0) ylab <- "y"

        dt <- .makedt(dt, annotations, ribosomal, reference)
        .bplot(dt)

    }
)

