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
                    ", or character vectors of KEGG/eggNOG annotations",
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
                nam[nam == ""] <- paste("subset", (1:nsubs)[nam == ""],
                                        sep = ".")
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
#' Calculate MILC values for every sequence in the given
#' \code{codonTable} object.
#'
#' @param cTobject A \code{codonTable} object.
#' @param subsets A (named) list of logical vectors, the length of each equal
#'    to \code{getlen(cTobject)}, i.e. the number of sequences in the set, or
#'    character vectors (of any length) containing KEGG/eggNOG annotations,
#'    or codonTable objects (of any length).
#' @param self Logical, if \code{TRUE} (default), CU statistic is also
#'    calculated against the average CU of the entire set of sequences.
#' @param ribosomal Logical, if \code{TRUE}, CU statistic is also calculated
#'    against the average CU of the ribosomal genes in the sequence set.
#' @inheritParams Biostrings::getGeneticCode
#' @param alt.init logical, whether to use alternative initiation codons.
#'    Default is \code{TRUE}.
#' @param stop.rm Logical, whether to remove stop codons. Default is
#'    \code{FALSE}.
#' @param filtering Character vector, one of \code{c("none", "soft", "hard")}.
#'    Specifies whether sequences shorther than some threshold value of
#'    length (in codons), \code{len.threshold}, should be excluded from
#'    calculations. If \code{"none"} (default), length of sequences is not
#'    checked, if \code{"soft"}, a warrning is printed if there are shorter
#'    sequences, and if \code{"hard"}, these sequences are excluded from
#'    calculation.
#' @param len.threshold Optional numeric, specifying sequence length,
#'    in codons, used for filtering.
#'
#'
#' @return A matrix with MILC values for every specified subset
#'    (\code{subsets}, \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of MILC, see
#'    \href{https://bit.ly/2GkT7qe}{Supek & Vlahovicek (2005)}.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#'
#' # calculate MILC distance to the average CU of the example DNA sequences
#' milc <- MILC(cT)
#' head(milc)
#'
#' # also calculate MILC distance to the average CU
#' # of ribosomal genes among the example DNA sequences
#' milc <- MILC(cT, ribosomal = TRUE)
#' head(milc)
#'
#' # calculate MILC distance to the average CU
#' # of the first 20 example DNA sequences
#' # (i.e. the first half of the example DNA set)
#' milc <- MILC(cT, self = FALSE,
#'              subsets = list(half = c(rep(TRUE, 20), rep(FALSE, 20))))
#'
#' # alternatively, you can specify codonTable as a subset
#' halfcT <- codonTable(codonCounts(cT)[1:20,])
#' milc2 <- MILC(cT, self = FALSE, subsets = list(half = halfcT))
#' all.equal(milc, milc2) # TRUE
#'
#' # filtering
#' MILC(cT, filtering = "hard", len.threshold = 80) # MILC for 9 sequences
#' sum(getlen(cT) > 80) # 9 sequences are longer than 80 codons
#' milc1 <- MILC(cT, filtering = "none") # no filtering
#' milc2 <- MILC(cT, filtering = "soft") # warning
#' all.equal(milc1, milc2) # TRUE
#'
#' # options for genetic code
#' milc <- MILC(cT, stop.rm = TRUE) # don't use stop codons in calculation
#' milc <- MILC(cT, alt.init = FALSE) # don't use alternative start codons
#' milc <- MILC(cT, id_or_name2 = "2") # use different genetic code, for help
#'                                     # see `?Biostrings::GENETIC_CODE`
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
            rowSums(2 * cTobject@counts * log(t(t(fc) / gc[, x])),
                    na.rm = TRUE) / len - cor
        })

    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### B, Karlin et al. 2001.
###
#' Calculate codon usage bias (B)
#'
#' Calculate B values for every sequence in the given
#' \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with B values for every specified subset (\code{subsets},
#'    \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of codon usage bias (B), see
#'    \href{https://bit.ly/2DRGdeb}{Karlin et al. (2001)}.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#'
#' # calculate B distance to the average CU of the example DNA sequences
#' b <- B(cT)
#' head(b)
#'
#' # also calculate B distance to the average CU
#' # of ribosomal genes among the example DNA sequences
#' b <- B(cT, ribosomal = TRUE)
#' head(b)
#'
#' # calculate B distance to the average CU
#' # of the first 20 example DNA sequences
#' # (i.e. the first half of the example DNA set)
#' b <- B(cT, self = FALSE,
#'        subsets = list(half = c(rep(TRUE, 20), rep(FALSE, 20))))
#'
#' # alternatively, you can specify codonTable as a subset
#' halfcT <- codonTable(codonCounts(cT)[1:20,])
#' b2 <- B(cT, self = FALSE, subsets = list(half = halfcT))
#' all.equal(b, b2) # TRUE
#'
#' # filtering
#' B(cT, filtering = "hard", len.threshold = 80) # B values for 9 sequences
#' sum(getlen(cT) > 80) # 9 sequences are longer than 80 codons
#' b1 <- B(cT, filtering = "none") # no filtering
#' b2 <- B(cT, filtering = "soft") # warning
#' all.equal(b1, b2) # TRUE
#'
#' # options for genetic code
#' b <- B(cT, stop.rm = TRUE) # don't use stop codons in calculation
#' b <- B(cT, alt.init = FALSE) # don't use alternative start codons
#' b <- B(cT, id_or_name2 = "2") # use different genetic code, for help
#'                               # see `?Biostrings::GENETIC_CODE`
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
#' Calculate MCB values for every sequence in the given
#' \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with MCB values for every specified subset
#'    (\code{subsets}, \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of maximum-likelihood codon bias (MCB), see
#'    \href{https://bit.ly/2GlMRyy}{Urrutia and Hurst (2001)}.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#'
#' # calculate MCB distance to the average CU of the example DNA sequences
#' mcb <- MCB(cT)
#' head(mcb)
#'
#' # also calculate MCB distance to the average CU
#' # of ribosomal genes among the example DNA sequences
#' mcb <- MCB(cT, ribosomal = TRUE)
#' head(mcb)
#'
#' # calculate MCB distance to the average CU
#' # of the first 20 example DNA sequences
#' # (i.e. the first half of the example DNA set)
#' mcb <- MCB(cT, self = FALSE,
#'            subsets = list(half = c(rep(TRUE, 20), rep(FALSE, 20))))
#'
#' # alternatively, you can specify codonTable as a subset
#' halfcT <- codonTable(codonCounts(cT)[1:20,])
#' mcb2 <- MCB(cT, self = FALSE, subsets = list(half = halfcT))
#' all.equal(mcb, mcb2) # TRUE
#'
#' # filtering
#' MCB(cT, filtering = "hard", len.threshold = 80) # MCB for 9 sequences
#' sum(getlen(cT) > 80) # 9 sequences are longer than 80 codons
#' mcb1 <- MCB(cT, filtering = "none") # no filtering
#' mcb2 <- MCB(cT, filtering = "soft") # warning
#' all.equal(mcb1, mcb2) # TRUE
#'
#' # options for genetic code
#' mcb <- MCB(cT, stop.rm = TRUE) # don't use stop codons in calculation
#' mcb <- MCB(cT, alt.init = FALSE) # don't use alternative start codons
#' mcb <- MCB(cT, id_or_name2 = "2") # use different genetic code, for help
#'                                   # see `?Biostrings::GENETIC_CODE`
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
            dt[which(is.na(dt), arr.ind = TRUE)] <- 0 # gc > 0
            dt[cTobject@counts<=0] <- 0 # counts > 0
            ba <- sapply(gCobject@cl, function(y)
                dt[, Reduce('+', .SD), .SDcols = y])
            rowSums(ba[, deg>1] * log10(csums[, deg>1]),
                    na.rm = TRUE) / rowSums(A)
        })
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ENC' Npvembre 2002.
###
#' Calculate effective number of codons prime (ENC').
#'
#' Calculate ENC' values for every sequence in the given
#' \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with ENC' values for every specified subset
#'    (\code{subsets}, \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of effective number of codons prime (ENC'), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/12140252}{Novembre (2002)}.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#'
#' # calculate ENC' distance to the average CU of the example DNA sequences
#' encp <- ENCprime(cT)
#' head(encp)
#'
#' # also calculate ENC' distance to the average CU
#' # of ribosomal genes among the example DNA sequences
#' encp <- ENCprime(cT, ribosomal = TRUE)
#' head(encp)
#'
#' # calculate ENC' distance to the average CU
#' # of the first 20 example DNA sequences
#' # (i.e. the first half of the example DNA set)
#' encp <- ENCprime(cT, self = FALSE,
#'                  subsets = list(half = c(rep(TRUE, 20), rep(FALSE, 20))))
#'
#' # alternatively, you can specify codonTable as a subset
#' halfcT <- codonTable(codonCounts(cT)[1:20,])
#' encp2 <- ENCprime(cT, self = FALSE, subsets = list(half = halfcT))
#' all.equal(encp, encp2) # TRUE
#'
#' # filtering
#' ENCprime(cT, filtering = "hard", len.threshold = 80) # ENC' for 9 sequences
#' sum(getlen(cT) > 80) # 9 sequences are longer than 80 codons
#' encp1 <- ENCprime(cT, filtering = "none") # no filtering
#' encp2 <- ENCprime(cT, filtering = "soft") # warning
#' all.equal(encp1, encp2) # TRUE
#'
#' # options for genetic code
#' encp <- ENCprime(cT, stop.rm = TRUE) # don't use stop codons in calculation
#' encp <- ENCprime(cT, alt.init = FALSE) # don't use alternative start codons
#' encp <- ENCprime(cT, id_or_name2 = "2") # use different genetic code,
#'                                         # see `?Biostrings::GENETIC_CODE`
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
            # if (any(gc==0, na.rm = TRUE)) {
            #     infs <- unique(which(dt == Inf, arr.ind = T)[,2])
            #     dt[, (infs) := 0]
            # }
            dt[which(is.na(dt), arr.ind = TRUE)] <- 0
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
#' Calculate ENC values for every sequence in the given \code{codonTable}
#' object.
#'
#' @inheritParams MILC
#'
#' @return A numeric vector with ENC values.
#'
#'    For definition of effective number of codons (ENC), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/2110097}{Wright (1990)}.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#'
#' # calculate ENC
#' enc <- ENC(cT)
#' head(enc)
#'
#' # filtering
#' ENC(cT, filtering = "hard", len.threshold = 80) # ENC for 9 sequences
#' sum(getlen(cT) > 80) # 9 sequences are longer than 80 codons
#' enc1 <- ENC(cT, filtering = "none") # no filtering
#' enc2 <- ENC(cT, filtering = "soft") # warning
#' all.equal(enc1, enc2) # TRUE
#'
#' # options for genetic code
#' enc <- ENC(cT, stop.rm = TRUE) # don't use stop codons in calculation
#' enc <- ENC(cT, alt.init = FALSE) # don't use alternative start codons
#' enc <- ENC(cT, id_or_name2 = "2") # use different genetic code, for help
#'                                   # see `?Biostrings::GENETIC_CODE`
#'
#' @rdname ENC
#' @export
setGeneric(
    name = "ENC",
    def = function(cTobject, id_or_name2 = "1",
                   alt.init = TRUE, stop.rm = TRUE,
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
        csums[which(csums <= 1, arr.ind = TRUE)] <- NA
        fa <- as.data.table((csums * pi - 1) / (csums - 1))
        .effNc(fa, gCobject@deg)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SCUO, Wan et al. 2004.
###
#' Calculate eorderliness of synonymous codon usage  (SCUO).
#'
#' Calculate SCUO values for every sequence in the given \code{codonTable}
#' object.
#'
#' @inheritParams MILC
#'
#' @return A numeric vector with SCUO values.
#'
#'    For definition of synonymous codon usage eorderliness (SCUO), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/15222899}{Wan et al. (2004)}.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#'
#' # calculate SCUO
#' scuo <- SCUO(cT)
#' head(scuo)
#'
#' # filtering
#' SCUO(cT, filtering = "hard", len.threshold = 80) # SCUO for 9 sequences
#' sum(getlen(cT) > 80) # 9 sequences are longer than 80 codons
#' scuo1 <- SCUO(cT, filtering = "none") # no filtering
#' scuo2 <- SCUO(cT, filtering = "soft") # warning
#' all.equal(scuo1, scuo2) # TRUE
#'
#' # options for genetic code
#' scuo <- SCUO(cT, stop.rm = TRUE) # don't use stop codons in calculation
#' scuo <- SCUO(cT, alt.init = FALSE) # don't use alternative start codons
#' scuo <- SCUO(cT, id_or_name2 = "2") # use different genetic code, for help
#'                                     # see `?Biostrings::GENETIC_CODE`
#'
#' @rdname SCUO
#' @export
setGeneric(
    name = "SCUO",
    def = function(cTobject, id_or_name2 = "1",
                   alt.init = TRUE, stop.rm = FALSE,
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
