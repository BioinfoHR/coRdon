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
            stop(paste0("Length of annotations vector, ", la,
                        ", differs from the number of sequences, ", nrow(dt),
                        "."))
        else dt[, annot := annotations]
    }

    dt[, genes := "other"]

    if (ribosomal){
        if (la == 0)
            stop("Ribosomal is TRUE, but no annotations were given!")
        else dt[annot %in% RPKOs, genes := "ribosomal"]
    }

    if (length(reference) != 0){
        if (is.null(names(reference)))
            stop("Reference must be a named list of length 1!")
        if (all(is.logical(unlist(reference)))) {
            r <- unlist(reference)
            n <- names(reference)
            dt[r, genes := n]
        } else if (all(is.character(unlist(reference)))){
            r <- unlist(reference)
            n <- names(reference)
            dt[annot %in% r, genes := n]
        }
    }
    dt
}
.bplot <- function(dt, size, alpha, guide){
        gp <- ggplot(dt, aes(x, y, colour = genes)) +
                geom_point(
                shape = 20, size = size, alpha = alpha, 
                show.legend = I(guide)
            ) + 
            guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
            theme_light()
        if (any(dt$genes != "other"))
        {
            alpha2 <- alpha + 0.3
            if (alpha2 > 1) alpha2 <- 1
            gp <- gp + geom_point(
                data = dt[genes != "other", ],
                alpha = alpha2, shape = 20, size = 1.5*size)
        }
        gp
}
#' Karlin B plot
#'
#' Plot distances of each gene's CU frequency to specified gene (sub)sets
#' (given by \code{x} and \code{y}).
#'
#' @param x,y Character, both must be in \code{colnames(data)}, or numeric
#'   vectors of CU statistic values for two subsets of genes. If numeric,
#'   the vectors must be of the same length.
#' @param data A matrix with CU statistic values for subsets of genes
#'   in columns.
#' @param annotations A character vector giving KO annotations for sequences
#'   for which the CU values were calculated, must be of length
#'   \code{nrow(data)}.
#' @param ribosomal Logical, whether to indicate ribosomal genes in the plot.
#'   Default is \code{FALSE}, if set to \code{TRUE}, then \code{annotation}
#'   must be given.
#' @param reference A named list of length 1, containing either a logical
#'   vector of \code{nrow(data)} of reference genes to be indicated
#'   on the plot, or a character vector (of any length) of the reference
#'   genes' anotations. If latter is the case, then \code{annotation}
#'   must be given.
#' @param size Numeric, indicating points' size
#' @param alpha Numeric, between 0 and 1, indicating points' transparency
#'   (default is 0.1).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' require(ggplot2)
#'
#' # calculate MILC distance to the average CU of the example DNA sequences,
#' # and to the average CU of ribosomal genes among the example DNA sequences
#' milc <- MILC(LD94, self = TRUE, ribosomal = TRUE)
#'
#' Bplot(x = "ribosomal", y = "self", data = milc,
#'       ribosomal = TRUE, annotations = getKO(LD94),
#'       size = 3) +
#'     labs(x = "MILC distance to ribosomal genes",
#'          y = "MILC distance to genes' average CU")
#'
#' @rdname Bplot
#' @export
setGeneric(
    name = "Bplot",
    def = function(
        x, y, data, 
        annotations = character(), ribosomal = FALSE, reference = list(),
        size = 1, alpha = 0.5
    )
    {
        standardGeneric("Bplot")
    }
)

#' @rdname Bplot
setMethod(
    f = "Bplot",
    signature = c(x = "character", y = "character", data = "matrix"),
    definition = function(
        x, y, data,
        annotations, ribosomal, reference, 
        size, alpha
    )
    {

        if (length(reference) == 0 &
            ribosomal == FALSE) guide <- FALSE
        else guide <- TRUE

        dt <- as.data.table(data[, c(x,y)])

        setnames(dt, c(x, y), c("x","y"))
        dt <- .makedt(dt, annotations, ribosomal, reference)
        .bplot(dt, size, alpha, guide)
    }
)
#' @rdname Bplot
setMethod(
    f = "Bplot",
    signature = c(x = "numeric", y = "numeric", data = "missing"),
    definition = function(
        x, y, data,
        annotations, ribosomal, reference, 
        size, alpha
    )
    {
        if (length(reference) == 0 &
            ribosomal == FALSE) guide <- FALSE
        else guide <- TRUE

        dt <- data.table(cbind(x, y))

        dt <- .makedt(dt, annotations, ribosomal, reference)
        .bplot(dt, size, alpha, guide)

    }
)

#' Intra-samples Karlin B plot
#'
#' Plot CU frequency distances between two samples
#' (given by \code{x} and \code{y}).
#'
#' @param x,y Objects of \code{codonTable} class.
#' @param names Character vector of length 2, giving names for samples.
#' @param variable A character, name of the function that will be used
#'   to calculate CU statistic values for plotting. Must be one of
#'   the following: \code{c("MILC", "B", "MCB", "ENCprime")}.
#' @param ribosomal Logical, whether to indicate ribosomal genes in the plot.
#'   Default is \code{FALSE}.
#' @param size Numeric, indicating points' size
#' @param alpha Numeric, between 0 and 1, indicating points' transparency
#'   (default is 0.1).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' require(ggplot2)
#' # calculate MILC distance to the average CU of the example DNA sequences,
#' # and to the average CU of ribosomal genes among the example DNA sequences
#' milc <- MILC(LD94, self = TRUE, ribosomal = TRUE)
#'
#' intraBplot(x = HD59, y = LD94, names = c("HD59", "LD94"),
#'            variable = "MILC", size = 3)
#'
#' @rdname intraBplot
#' @export
setGeneric(
    name = "intraBplot",
    def = function(
        x, y, names = c("x", "y"), variable, ribosomal = FALSE, 
        size = 1, alpha = 0.5
    )
    {
        standardGeneric("intraBplot")
    }
)

#' @rdname intraBplot
setMethod(
    f = "intraBplot",
    signature = c(x = "codonTable", y = "codonTable"),
    definition = function(x, y, names, variable, ribosomal, size, alpha)
    {
        sl <- list(x, y)
        names(sl) <- names

        if (variable == "MILC") {
            bpf <- function(...) MILC(...)
        } else if (variable == "B") {
            bpf <- function(...) B(...)
        } else if (variable == "MCB") {
            bpf <- function(...) MCB(...)
        } else if (variable == "ENCprime") {
            bpf <- function(...) ENCprime(...)
        }
        dtx <- as.data.table(bpf(x, self = FALSE, subsets = sl))
        dty <- as.data.table(bpf(y, self = FALSE, subsets = sl))
        rsl <- list(dtx, dty)
        names(rsl) <- names
        dt <- rbindlist(rsl, idcol = "sample")

        gp <- ggplot(dt, aes(get(names[1]), get(names[2]), colour = sample)) +
            geom_point(shape = 20, size = size, alpha = alpha) +
            guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
            labs(x = names[1], y = names[2]) +
            theme_light()

        if (ribosomal) {
            alpha2 <- alpha + 0.3
            if (alpha2 > 1) alpha2 <- 1
            rows <- c(getKO(x) %in% RPKOs, getKO(y) %in% RPKOs)
            gp <- gp + geom_point(
                data = dt[rows, ], 
                aes(get(names[1]), get(names[2]), colour = sample),
                alpha = alpha2, shape = 20, size = 1.5*size)
        }
        gp

    }
)
