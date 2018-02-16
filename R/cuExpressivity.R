#' @include codonTable-class.R
#' @include genCode-class.R
#' @include functions.R
#' @include codonUsage.R
#' @import data.table

setGeneric(
    name = "MELP",
    def = function(cTobject, subsets = list(), ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE) {
        standardGeneric("MELP")
    }
)

setMethod(
    f = "MELP",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, subsets, ribosomal,
                          id_or_name2, alt.init, stop.rm) {

        milcs <- MILC(cTobject, subsets, self = TRUE, ribosomal,
             id_or_name2, alt.init, stop.rm)
        sapply(colnames(milcs)[-1], function(y) milcs[,"self"] / milcs[, y])
    }
)

setGeneric(
    name = "E",
    def = function(cTobject, subsets = list(), ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE) {
        standardGeneric("E")
    }
)

setMethod(
    f = "E",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, subsets, ribosomal,
                          id_or_name2, alt.init, stop.rm) {

        b <- B(cTobject, subsets, self = TRUE, ribosomal,
               id_or_name2, alt.init, stop.rm)
        sapply(colnames(b)[-1], function(y) b[,"self"] / b[, y])
    }
)

setGeneric(
    name = "CAI",
    def = function(cTobject, subsets = list(), ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE) {
        standardGeneric("CAI")
    }
)

setMethod(
    f = "CAI",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, subsets, ribosomal,
                          id_or_name2, alt.init, stop.rm) {

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

setGeneric(
    name = "Fop",
    def = function(cTobject, subsets = list(), ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE) {
        standardGeneric("Fop")
    }
)

setMethod(
    f = "Fop",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, subsets, ribosomal,
                          id_or_name2, alt.init, stop.rm)  {

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

setGeneric(
    name = "GCB",
    def = function(cTobject, seed, ribosomal = FALSE, perc = 0.05,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE) {
        standardGeneric("GCB")
    }
)

setMethod(
    f = "GCB",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, seed, ribosomal, perc,
                          id_or_name2, alt.init, stop.rm)  {

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
