#' @include codonTable-class.R
#' @include genCode-class.R
#' @include functions.R
#' @import data.table

setGeneric(
    name = "expectedCU",
    def = function(cTobject, gCobject, subsets, self, ribosomal){
        standardGeneric("expectedCU")
    }
)

setMethod(
    f = "expectedCU",
    signature = c("codonTable", "genCode", "list", "logical", "logical"),
    definition = function(cTobject, gCobject, subsets, self, ribosomal) {
        ns <- length(cTobject@len)
        if (!is.list(subsets))
            stop(
                paste0(
                    "Subsets must be either a (named) list of logical vectors",
                    ", each of length ",
                    ns,
                    ", or character vectors containing KEGG/eggNOG annotations",
                    ", or codonTable objects (of any length)."
                )
            )
        if (length(subsets) != 0) {
            ok <- sapply(subsets, function(x) {
                all(is.vector(x, mode = "logical"), length(x) == ns) |
                    is.vector(x, mode = "character") |
                    all(inherits(x, "codonTable"), nrow(x) > 0)
            })
            stopifnot(ok)
            nam <- names(subsets)
            nsubs <- length(subsets)
            if (is.null(nam)) {
                nam <- paste("subset", 1:nsubs, sep = ".")
            } else {
                nam[nam == ""] <- paste("subset", (1:nsubs)[nam == ""], sep = ".")
            }
            names(subsets) <- make.names(nam, unique = TRUE)
        }
        if (ribosomal) {
            ribKOs <- cTobject@KO %in% RPKOs
            if (any(ribKOs))
                subsets <- c(list(ribosomal = ribKOs), subsets)
        }
        if (self) {
            self_set <- rep.int(TRUE, ns)
            subsets <- c(list(self = c(self_set)), subsets)
        }
        sapply(subsets, function(s) {
            if (all(class(s) %in% c("logical", "character")))
                .normSetFrequencies(subset(cTobject, s), gCobject)
            else if (all(class(s) == "codonTable"))
                .normSetFrequencies(s, gCobject)
        })
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MILC, Supek and Vlahovicek 2005.
###
#' @export
setGeneric(
    name = "MILC",
    def = function(cTobject, subsets = list(), self = TRUE, ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("MILC")
    }
)
#' @export
setMethod(
    f = "MILC",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, subsets, self, ribosomal,
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

        len <- cTobject@len
        csums <- .countsbyaa(cTobject, gCobject)
        fc <- .normFrequencies(cTobject, gCobject)
        gc <- expectedCU(cTobject, gCobject, subsets, self, ribosomal)
        cor <- as.numeric(((csums > 0) %*% (gCobject@deg - 1)) / len - 0.5)
        sapply(colnames(gc), function(x) {
            rowSums(2 * cTobject@counts * log(t(t(fc) / gc[, x])), na.rm = TRUE) /
                len - cor
        })

    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### B, Karlin and Mrazek 1998.
###
#' @export
setGeneric(
    name = "B",
    def = function(cTobject, subsets = list(), self = TRUE, ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("B")
    }
)
#' @export
setMethod(
    f = "B",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, subsets, self, ribosomal,
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

        pa <- .countsbyaa(cTobject, gCobject) / cTobject@len
        fc <- .normFrequencies(cTobject, gCobject)
        gc <- expectedCU(cTobject, gCobject, subsets, self, ribosomal)
        sapply(colnames(gc), function(x) {
            dt <- as.data.table(t(abs(t(fc) - gc[, x])))
            ba <- sapply(gCobject@cl, function(x)
                dt[, Reduce('+', .SD), .SDcols = x])
            rowSums(pa * ba, na.rm = TRUE)
        })
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MCB, Urutia and Hurst 2001.
###
#' @export
setGeneric(
    name = "MCB",
    def = function(cTobject, subsets = list(), self = TRUE, ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("MCB")
    }
)
#' @export
setMethod(
    f = "MCB",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, subsets, self, ribosomal,
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
        csums <- .countsbyaa(cTobject, gCobject)
        A <- csums[, deg>1] > 0
        fc <- .normFrequencies(cTobject, gCobject)
        gc <- expectedCU(cTobject, gCobject, subsets, self, ribosomal)
        sapply(colnames(gc), function(x) {
            dt <- as.data.table(t((t(fc) - gc[, x])^2 / gc[, x]))
            dt[which(is.na(dt), arr.ind = T)] <- 0 # gc > 0
            dt[cTobject@counts<=0] <- 0 # counts > 0
            ba <- sapply(gCobject@cl, function(y)
                dt[, Reduce('+', .SD), .SDcols = y])
            rowSums(ba[, deg>1] * log10(csums[, deg>1]), na.rm = T) / rowSums(A)
        })
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ENC' Npvembre 2002.
###
#' @export
setGeneric(
    name = "ENCprime",
    def = function(cTobject, subsets = list(), self = TRUE, ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = TRUE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("ENCprime")
    }
)
#' @export
setMethod(
    f = "ENCprime",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, subsets, self, ribosomal,
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
        csums <- .countsbyaa(cTobject, gCobject)
        fc <- .normFrequencies(cTobject, gCobject)
        gc <- expectedCU(cTobject, gCobject, subsets, self, ribosomal)
        sapply(colnames(gc), function(x) {
            # chi squared
            dt <- as.data.table(t((t(fc) - gc[, x])^2 / gc[, x]))
            # if (any(gc==0, na.rm = T)) {
            #     infs <- unique(which(dt == Inf, arr.ind = T)[,2])
            #     dt[, (infs) := 0]
            # }
            dt[which(is.na(dt), arr.ind = T)] <- 0
            chisum <- sapply(gCobject@cl, function(x)
                dt[, Reduce('+', .SD), .SDcols = x])
            chisq <- csums * chisum
            # homozygosity
            fa <- as.data.table(t((t(chisq + csums) - deg) /
                                      (deg * t(csums - 1))))
            fa[csums<5] <- NA
            .effNc(fa, deg)
        })
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ENC, Wright 1990.
###
#' @export
setGeneric(
    name = "ENC",
    def = function(cTobject, id_or_name2 = "1", alt.init = TRUE, stop.rm = TRUE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("ENC")
    }
)
#' @export
setMethod(
    f = "ENC",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, id_or_name2, alt.init, stop.rm,
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
        cl <- gCobject@cl
        counts <- as.data.table(cTobject@counts)
        csums <- .countsbyaa(cTobject, gCobject)
        freqs <- sapply(seq_along(cl), function(x)
            counts[,.SD/csums[,x], .SDcols = cl[[x]]])
        pi <- sapply(freqs, function(x) rowSums(x^2, na.rm = TRUE))
        # remove AA with count <= 1
        csums[which(csums <= 1, arr.ind = T)] <- NA
        fa <- as.data.table((csums * pi - 1) / (csums - 1))
        .effNc(fa, gCobject@deg)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SCUO, Wan et al. 2004.
###
#' @export
setGeneric(
    name = "SCUO",
    def = function(cTobject, id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("SCUO")
    }
)
#' @export
setMethod(
    f = "SCUO",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject, id_or_name2, alt.init, stop.rm,
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

        cl <- gCobject@cl
        counts <- as.data.table(cTobject@counts)
        csums <- .countsbyaa(cTobject, gCobject)
        freqs <- sapply(seq_along(cl), function(x)
            counts[,.SD/csums[,x], .SDcols = cl[[x]]])
        Ha <- sapply(freqs, function(x) rowSums(-x*log10(x), na.rm = TRUE))
        Hmax <- log10(sapply(cl, length))
        Oa <- t((Hmax - t(Ha)) / Hmax)
        Fa <- csums / rowSums(csums[,gCobject@deg>1])
        rowSums(Oa * Fa, na.rm = TRUE)
    }
)
