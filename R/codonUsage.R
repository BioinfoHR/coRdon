#' @include codonTable-class.R
#' @include genCode-class.R
#' @include functions.R
#' @import data.table
NULL

setGeneric(
    name = "expectedCU",
    def = function(cTobject, gCobject, subsets, self, ribosomal){
        standardGeneric("expectedCU")
    }
)

setMethod(
    f = "expectedCU",
    signature = c(cTobject = "codonTable", subsets = "list"),
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
#' Calculate codon usage Measure Independent of Length and Composition (MILC)
#'
#' Calculate MILC values for every sequence in the given \code{codonTable} object.
#'
#' @param cTobject A \code{codonTable} object.
#' @param subsets A (named) list of logical vectors, the length of each equal
#'    to \code{getlen(cTobject)}, i.e. the number of sequences in the set, or
#'    character vectors (of any length) containing KEGG/eggNOG annotations,
#'    or codonTable objects (of any length).
#' @param self Logical, if \code{TRUE} (default), CU statistic is also calculated
#'    against the average CU of the entire set of sequences.
#' @param ribosomal Logical, if \code{TRUE}, CU statistic is also calculated
#'    against the average CU of the ribosomal genes in the sequence set.
#' @inheritParams Biostrings::getGeneticCode
#' @param alt.init logical, whether to use alternative initiation codons.
#'    Default is \code{TRUE}.
#' @param stop.rm Logical, whether to remove stop codons. Default is \code{FALSE}.
#' @param filtering Character vector, one of \code{c("none", "soft", "hard")}.
#'    Specifies whether sequences shorther than some threshold value of length
#'    (in codons), \code{len.threshold}, should be excluded from calculations. If
#'    \code{"none"} (default), length of sequences is not checked, if \code{"soft"},
#'    a warrning is printed if there are shorter sequences, and if \code{"hard"},
#'    these sequences are excluded.
#' @param len.threshold Optional numeric, specifying sequence length, in codons,
#'    used for filtering.
#'
#'
#' @return A matrix with MILC values for every specified subset (\code{subsets},
#'    \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of MILC, see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/16029499}{Supek & Vlahovicek (2005)}.
#'
#' @rdname MILC
#' @export
setGeneric(
    name = "MILC",
    def = function(cTobject, subsets = list(), self = TRUE, ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("MILC")
    }
)

#' @rdname MILC
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
### B, Karlin et al. 2001.
###
#' Calculate codon usage bias (B)
#'
#' Calculate B values for every sequence in the given \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with B values for every specified subset (\code{subsets},
#'    \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of codon usage bias (B), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/11489855}{Karlin et al. (2001)}.
#'
#' @rdname B
#' @export
setGeneric(
    name = "B",
    def = function(cTobject, subsets = list(), self = TRUE, ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("B")
    }
)

#' @rdname B
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
#' Calculate maximum-likelihood codon bias (MCB)
#'
#' Calculate MCB values for every sequence in the given \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with MCB values for every specified subset (\code{subsets},
#'    \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of maximum-likelihood codon bias (MCB), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/2110097}{Urrutia and Hurst (2001)}.
#'
#' @rdname MCB
#' @export
setGeneric(
    name = "MCB",
    def = function(cTobject, subsets = list(), self = TRUE, ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("MCB")
    }
)

#' @rdname MCB
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
#' Calculate effective number of codons prime (ENC').
#'
#' Calculate ENC' values for every sequence in the given \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with ENC' values for every specified subset (\code{subsets},
#'    \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of effective number of codons prime (ENC'), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/12140252}{Novembre (2002)}.
#'
#' @rdname ENCprime
#' @export
setGeneric(
    name = "ENCprime",
    def = function(cTobject, subsets = list(), self = TRUE, ribosomal = FALSE,
                   id_or_name2 = "1", alt.init = TRUE, stop.rm = TRUE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("ENCprime")
    }
)

#' @rdname ENCprime
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
#' Calculate effective number of codons (ENC).
#'
#' Calculate ENC values for every sequence in the given \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A numeric vector with ENC values.
#'
#'    For definition of effective number of codons (ENC), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/2110097}{Wright (1990)}.
#'
#' @rdname ENC
#' @export
setGeneric(
    name = "ENC",
    def = function(cTobject, id_or_name2 = "1", alt.init = TRUE, stop.rm = TRUE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("ENC")
    }
)

#' @rdname ENC
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
#' Calculate eorderliness of synonymous codon usage  (SCUO).
#'
#' Calculate SCUO values for every sequence in the given \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A numeric vector with SCUO values.
#'
#'    For definition of synonymous codon usage eorderliness (SCUO), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/15222899}{Wan et al. (2004)}.
#'
#' @rdname SCUO
#' @export
setGeneric(
    name = "SCUO",
    def = function(cTobject, id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE,
                   filtering = "none", len.threshold = 80) {
        standardGeneric("SCUO")
    }
)

#' @rdname SCUO
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
