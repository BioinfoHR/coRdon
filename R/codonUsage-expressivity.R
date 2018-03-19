#' @include codonTable-class.R
#' @include genCode-class.R
#' @include functions.R
#' @include codonUsage.R
#' @import data.table
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MELP, Supek and Vlahovicek 2005.
###
#' @export
setGeneric(
    name = "MELP",
    def = function(cTobject, subsets = list(), ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("MELP")
    }
)
#' Calculate MILC-based Expression Level Predictor (MELP).
#'
#' Calculate MELP values for every sequence in the given \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with MELP values for every specified subset (\code{subsets},
#'    \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of MELP, see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/16029499}{Supek & Vlahovicek (2005)}.
#'
#' @name MELP
#' @export
setMethod(
    f = "MELP",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, subsets, ribosomal,
                          id_or_name2, alt.init, stop.rm,
                          filtering, len.threshold) {

        if (filtering == "hard") {
            cTobject <- subset(cTobject, cTobject@len > len.threshold)
        } else if (filtering == "soft") {
            if (any(cTobject@len < len.threshold))
                warning("Some sequences have below-threshold length!")
        } else if (filtering == "none") NULL

        milcs <- MILC(cTobject, subsets, self = TRUE, ribosomal,
             id_or_name2, alt.init, stop.rm)
        sapply(colnames(milcs)[-1], function(y) milcs[,"self"] / milcs[, y])
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### E, Karlin and Mrazek 2000
###
#' @export
setGeneric(
    name = "E",
    def = function(cTobject, subsets = list(), ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("E")
    }
)
#' Calculate gene expression measure (E).
#'
#' Calculate E values for every sequence in the given \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with E values for every specified subset (\code{subsets},
#'    \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of gene expression measure (E), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/10960111}{Karlin and Mrazek (2000)}.
#'
#' @name E
#' @export
setMethod(
    f = "E",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, subsets, ribosomal,
                          id_or_name2, alt.init, stop.rm,
                          filtering, len.threshold) {

        if (filtering == "hard") {
            cTobject <- subset(cTobject, cTobject@len > len.threshold)
        } else if (filtering == "soft") {
            if (any(cTobject@len < len.threshold))
                warning("Some sequences have below-threshold length!")
        } else if (filtering == "none") NULL

        b <- B(cTobject, subsets, self = TRUE, ribosomal,
               id_or_name2, alt.init, stop.rm)
        sapply(colnames(b)[-1], function(y) b[,"self"] / b[, y])
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### CAI, Sharp & Lee 1987
###
#' @export
setGeneric(
    name = "CAI",
    def = function(cTobject, subsets = list(), ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("CAI")
    }
)
#' Calculate Codon Adaptation Index (CAI).
#'
#' Calculate CAI values for every sequence in the given \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with CAI values for every specified subset (\code{subsets},
#'    \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of Codon Adaptation Index (CAI), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/3547335}{Sharp and Li (1987)}.
#'
#' @name CAI
#' @export
setMethod(
    f = "CAI",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, subsets, ribosomal,
                          id_or_name2, alt.init, stop.rm,
                          filtering, len.threshold) {

        if (filtering == "hard") {
            cTobject <- subset(cTobject, cTobject@len > len.threshold)
        } else if (filtering == "soft") {
            if (any(cTobject@len < len.threshold))
                warning("Some sequences have below-threshold length!")
        } else if (filtering == "none") NULL

        gCobject <- genCode(id_or_name2, alt.init, stop.rm)
        if (stop.rm) {
            cTobject@counts <- cTobject@counts[,gCobject@nostops]
            cTobject@len <- rowSums(cTobject@counts[,gCobject@nostops])
        }
        deg <- gCobject@deg
        aa <- unlist(gCobject@ctab[, "AA"], use.names = FALSE)
        counts <- as.data.table(cTobject@counts)
        csums <- .countsbyaa(cTobject, gCobject)
        fc <- .normFrequencies(cTobject, gCobject)
        gc <- expectedCU(cTobject, gCobject, subsets, self = FALSE, ribosomal)
        gcmax <-
            apply(gc, 2, function(x) # freq of the most frequent codons by aa
                by(x, droplevels(aa), max))
        gcmax[gcmax == 0] <- 0.5 # account for codons that are not used in ref. set
        x <- as.integer(droplevels(aa))
        gcm <- sapply(colnames(gc), function(y) # freq of the most frequent codons
            gc[, y] <- gcmax[x, y])
        nondeg <- match(which(deg == 1), as.integer(aa))
        counts[, (nondeg) := NA] # remove codons for non-degenerate aa
        sapply(colnames(gcm), function(y)
            exp(rowSums(counts * t(log(t(fc) / gcm[, y])), na.rm = TRUE) /
                    rowSums(counts, na.rm = T)))

    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Fop, Ikemura 1981
###
#' @export
setGeneric(
    name = "Fop",
    def = function(cTobject, subsets = list(), ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("Fop")
    }
)
#' Calculate frequency of optimal codons (Fop).
#'
#' Calculate Fop values for every sequence in the given \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with Fop values for every specified subset (\code{subsets},
#'    \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of frequency of optimal codons (Fop), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/6175758}{Ikemura (1981)}.
#'
#' @name Fop
#' @export
setMethod(
    f = "Fop",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, subsets, ribosomal,
                          id_or_name2, alt.init, stop.rm,
                          filtering, len.threshold)  {

        if (filtering == "hard") {
            cTobject <- subset(cTobject, cTobject@len > len.threshold)
        } else if (filtering == "soft") {
            if (any(cTobject@len < len.threshold))
                warning("Some sequences have below-threshold length!")
        } else if (filtering == "none") NULL

        gCobject <- genCode(id_or_name2, alt.init, stop.rm)
        if (stop.rm) {
            cTobject@counts <- cTobject@counts[,gCobject@nostops]
            cTobject@len <- rowSums(cTobject@counts[,gCobject@nostops])
        }
        deg <- gCobject@deg
        aa <- unlist(gCobject@ctab[, "AA"], use.names = FALSE)
        counts <- as.data.table(cTobject@counts)
        csums <- .countsbyaa(cTobject, gCobject)
        fc <- .normFrequencies(cTobject, gCobject)
        gc <- expectedCU(cTobject, gCobject, subsets, self = FALSE, ribosomal)
        gcmax <-
            apply(gc, 2, function(x) # freq of the most frequent codons by aa
                by(x, droplevels(aa), max))
        gcmax[gcmax == 0] <- 0.5 # account for codons that are not used in ref. set
        x <- as.integer(droplevels(aa))
        gcm <- sapply(colnames(gc), function(y) # freq of the most frequent codons
            gc[, y] <- gcmax[x, y])
        nondeg <- match(which(deg == 1), as.integer(aa))
        counts[, (nondeg) := NA] # remove codons for non-degenerate aa
        sapply(colnames(gcm), function(y){
            ra <- t(fc) / gcm[, y] # relative adaptivenes of codons
            cnt <- as.matrix(counts)
            cnt[which(t(ra) < 0.9, arr.ind = T)] <- NA # as implemented in INCA
            rowSums(cnt, na.rm = T) / rowSums(counts, na.rm = T)
        })
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GCB, Merkl 2003
###
#' @export
setGeneric(
    name = "GCB",
    def = function(cTobject, seed, ribosomal = FALSE, perc = 0.05,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("GCB")
    }
)
#' Calculate gene codon bias (GCB).
#'
#' Calculate GCB values for every sequence in the given \code{codonTable} object.
#'
#' @inheritParams MILC
#' @param seed A logical vector, of the length equal to \code{getlen(cTobject)},
#'    or a character vector (of any length) containing KEGG/eggNOG annotations,
#'    or a codonTable object (of any length). Indicates a set of genes, or their CU,
#'    to be used as a target in the first iteration of the algorithm.
#' @param ribosomal Logical, if \code{TRUE}, ribosomal genes are used as a seed.
#' @param perc percent of top ranking genes to be used as a target set for the
#'    next iteration
#'
#' @return A numeric vector with GCB values.
#'
#'    For definition of gene codon bias (GCB), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/6175758}{Ikemura (1981)}.
#'
#' @name GCB
#' @export
setMethod(
    f = "GCB",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, seed, ribosomal, perc,
                          id_or_name2, alt.init, stop.rm,
                          filtering, len.threshold)  {

        if (filtering == "hard") {
            cTobject <- subset(cTobject, cTobject@len > len.threshold)
        } else if (filtering == "soft") {
            if (any(cTobject@len < len.threshold))
                warning("Some sequences have below-threshold length!")
        } else if (filtering == "none") NULL

        if (length(seed) == 0 & ribosomal == FALSE)
            stop("Seed is not specified!")

        gCobject <- genCode(id_or_name2, alt.init, stop.rm)
        if (stop.rm) {
            cTobject@counts <- cTobject@counts[,gCobject@nostops]
            cTobject@len <- rowSums(cTobject@counts[,gCobject@nostops])
        }

        len <- cTobject@len
        counts <- as.data.table(cTobject@counts)
        fc <- .normFrequencies(cTobject, gCobject)
        gc <- expectedCU(cTobject, gCobject, list(seed=seed), self = FALSE, ribosomal)
        gcb_prev <- numeric(nrow(counts))
        iter <- 0
        repeat {
            cb <- log(gc / colMeans(fc, na.rm = T))
            cb[gc == 0] <- -5
            gcb <- rowSums(t(t(counts) * as.numeric(cb)), na.rm = T) / len
            diff <- all(gcb == gcb_prev)
            if (diff | iter > 6)
                break
            else {
                iter <- iter + 1
                gcb_prev <- gcb
                top <- order(gcb, decreasing = T)[1:(perc * length(gcb))]
                gc <- expectedCU(cTobject, gCobject,
                                 list(s = codonTable(counts[top,])), self = FALSE, ribosomal)
            }
        }
        return(gcb)
    }
)
