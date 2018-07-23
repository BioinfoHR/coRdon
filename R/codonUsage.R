#' @include codonTable-class.R
#' @include genCode-class.R
#' @include functions.R
#' @import data.table
NULL

setGeneric(
    name = "expectedCU",
    def = function(cTobject,
                    gCobject,
                    subsets,
                    self,
                    ribosomal)
    {
        standardGeneric("expectedCU")
    }
)

setMethod(
    f = "expectedCU",
    signature = c(cTobject = "codonTable", subsets = "list"),
    definition = function(cTobject,
                            gCobject,
                            subsets,
                            self,
                            ribosomal)
    {
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
            ok <- vapply(subsets, function(x) {
                all(is.vector(x, mode = "logical"), length(x) == ns) |
                    is.vector(x, mode = "character") |
                    all(inherits(x, "codonTable"), nrow(x) > 0)
            }, logical(length = 1))
            stopifnot(ok)
            nam <- names(subsets)
            nsubs <- length(subsets)
            if (is.null(nam)) {
                nam <- paste("subset", seq_len(nsubs), sep = ".")
            } else {
                nam[nam == ""] <- paste("subset", seq_len(nsubs)[nam == ""],
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
        vapply(subsets, function(s) {
            if (all(class(s) %in% c("logical", "character")))
                .normSetFrequencies(subset(cTobject, s), gCobject)
            else if (all(is(s, "codonTable")))
                .normSetFrequencies(s, gCobject)
        }, numeric(length = nrow(gCobject@ctab)))
    }
)

#' Calculate CU measures.
#'
#' Calculate values of the codon usage (CU) measure 
#' for every sequence in the given \code{codonTable} object.
#' The following methods are implemented:
#'  \code{MILC}, Measure Independent of Length and Composition  
#'  \href{https://bit.ly/2GkT7qe}{Supek & Vlahovicek (2005)}, 
#'  \code{B}, codon usage bias (B) 
#'  \href{https://bit.ly/2DRGdeb}{Karlin et al. (2001)}, 
#'  \code{ENC}, effective number of codons (ENC) 
#'  \href{https://www.ncbi.nlm.nih.gov/pubmed/2110097}{Wright (1990)}.
#'  \code{ENCprime}, effective number of codons prime (ENC') 
#'  \href{https://www.ncbi.nlm.nih.gov/pubmed/12140252}{Novembre (2002)}, 
#'  \code{MCB}, maximum-likelihood codon bias (MCB) 
#'  \href{https://bit.ly/2GlMRyy}{Urrutia and Hurst (2001)}, 
#'  \code{SCUO}, synonymous codon usage eorderliness (SCUO) 
#'  \href{https://www.ncbi.nlm.nih.gov/pubmed/15222899}{Wan et al. (2004)}.
#'  
#' @param cTobject A \code{codonTable} object.
#' @param subsets A (named) list of logical vectors, the length of each equal
#'    to \code{getlen(cTobject)}, i.e. the number of sequences in the set, or
#'    character vectors (of any length) containing KEGG/eggNOG annotations,
#'    or codonTable objects (of any length). 
#'    Not used  for \code{ENC}, \code{SCUO} and \code{GCB} calculations.
#' @param self Logical, if \code{TRUE} (default), CU statistic is also
#'    calculated against the average CU of the entire set of sequences.
#'    Not used for \code{ENC}, \code{SCUO} and \code{GCB} calculations.
#' @param ribosomal Logical, if \code{TRUE}, CU statistic is also calculated
#'    against the average CU of the ribosomal genes in the sequence set.
#'    Not used  for \code{ENC} and \code{SCUO} calculations.
#'    For GCB calculations, if \code{TRUE}, ribosomal genes are used 
#'    as a seed, and if \code{FALSE}  (default), \code{seed} 
#'    has to be specified.
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
#' @return A matrix or a numeric vector with CU measure values. 
#'    For \code{MILC}, \code{B}, \code{ENCprime}, the matrix has a column 
#'    with values for every specified subset
#'    (\code{subsets}, \code{self}, \code{ribosomal}).
#'    A numeric vector for \code{ENC} and \code{SCUO}.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#'
#' # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' # In the examples below, MILC values are calculated for all sequences; 
#' # B and ENCprime can be caluclated in the same way.
#' # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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
#' # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' # In the examples below, ENC values are calculated for all sequences; 
#' # SCUO values can be caluclated in the same way.
#' # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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
#' @name codonUsage
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MILC, Supek and Vlahovicek 2005.
###
#' @rdname codonUsage
#' @export
setGeneric(
    name = "MILC",
    def = function(cTobject,
                    subsets = list(),
                    self = TRUE,
                    ribosomal = FALSE,
                    id_or_name2 = "1",
                    alt.init = TRUE,
                    stop.rm = FALSE,
                    filtering = "none",
                    len.threshold = 80)
    {
        standardGeneric("MILC")
    }
)

#' @rdname codonUsage
setMethod(
    f = "MILC",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject,
                            subsets,
                            self,
                            ribosomal,
                            id_or_name2,
                            alt.init,
                            stop.rm,
                            filtering,
                            len.threshold) {
        if (filtering == "hard") {
            cTobject <- subset(cTobject, cTobject@len > len.threshold)
        } else if (filtering == "soft") {
            if (any(cTobject@len < len.threshold))
                warning("Some sequences have below-threshold length!")
        } else if (filtering == "none")
            NULL

        gCobject <- genCode(id_or_name2, alt.init, stop.rm)
        if (stop.rm) {
            cTobject@counts <- cTobject@counts[, gCobject@nostops]
            cTobject@len <-
                rowSums(cTobject@counts[, gCobject@nostops])
        }

        len <- cTobject@len
        csums <- .countsbyaa(cTobject, gCobject)
        fc <- .normFrequencies(cTobject, gCobject)
        gc <-
            expectedCU(cTobject, gCobject, subsets, self, ribosomal)
        cor <-
            as.numeric(((csums > 0) %*% (gCobject@deg - 1)) / len - 0.5)
        vapply(colnames(gc), function(x) {
            rowSums(2 * cTobject@counts * log(t(t(fc) / gc[, x])),
                    na.rm = TRUE) / len - cor
        }, numeric(length = nrow(cTobject@counts)))

    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### B, Karlin et al. 2001.
###
#' @rdname codonUsage
#' @export
setGeneric(
    name = "B",
    def = function(cTobject,
                    subsets = list(),
                    self = TRUE,
                    ribosomal = FALSE,
                    id_or_name2 = "1",
                    alt.init = TRUE,
                    stop.rm = FALSE,
                    filtering = "none",
                    len.threshold = 80)
    {
        standardGeneric("B")
    }
)

#' @rdname codonUsage
setMethod(
    f = "B",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject,
                            subsets,
                            self,
                            ribosomal,
                            id_or_name2,
                            alt.init,
                            stop.rm,
                            filtering,
                            len.threshold)
    {
        if (filtering == "hard") {
            cTobject <- subset(cTobject, cTobject@len > len.threshold)
        } else if (filtering == "soft") {
            if (any(cTobject@len < len.threshold))
                warning("Some sequences have below-threshold length!")
        } else if (filtering == "none")
            NULL

        gCobject <- genCode(id_or_name2, alt.init, stop.rm)
        if (stop.rm) {
            cTobject@counts <- cTobject@counts[, gCobject@nostops]
            cTobject@len <-
                rowSums(cTobject@counts[, gCobject@nostops])
        }

        pa <- .countsbyaa(cTobject, gCobject) / cTobject@len
        fc <- .normFrequencies(cTobject, gCobject)
        gc <-
            expectedCU(cTobject, gCobject, subsets, self, ribosomal)
        vapply(colnames(gc), function(x) {
            dt <- as.data.table(t(abs(t(fc) - gc[, x])))
            ba <- vapply(gCobject@cl, function(x)
                dt[, Reduce('+', .SD), .SDcols = x],
                numeric(length = nrow(dt)))
            rowSums(pa * ba, na.rm = TRUE)
        }, numeric(length = nrow(cTobject@counts)))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MCB, Urutia and Hurst 2001.
###
#' @rdname codonUsage
#' @export
setGeneric(
    name = "MCB",
    def = function(cTobject,
                    subsets = list(),
                    self = TRUE,
                    ribosomal = FALSE,
                    id_or_name2 = "1",
                    alt.init = TRUE,
                    stop.rm = FALSE,
                    filtering = "none",
                    len.threshold = 80)
    {
        standardGeneric("MCB")
    }
)

#' @rdname codonUsage
setMethod(
    f = "MCB",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject,
                            subsets,
                            self,
                            ribosomal,
                            id_or_name2,
                            alt.init,
                            stop.rm,
                            filtering,
                            len.threshold)
    {
        if (filtering == "hard") {
            cTobject <- subset(cTobject, cTobject@len > len.threshold)
        } else if (filtering == "soft") {
            if (any(cTobject@len < len.threshold))
                warning("Some sequences have below-threshold length!")
        } else if (filtering == "none")
            NULL

        gCobject <- genCode(id_or_name2, alt.init, stop.rm)
        if (stop.rm) {
            cTobject@counts <- cTobject@counts[, gCobject@nostops]
            cTobject@len <-
                rowSums(cTobject@counts[, gCobject@nostops])
        }

        deg <- gCobject@deg
        csums <- .countsbyaa(cTobject, gCobject)
        A <- csums[, deg > 1] > 0
        fc <- .normFrequencies(cTobject, gCobject)
        gc <-
            expectedCU(cTobject, gCobject, subsets, self, ribosomal)
        vapply(colnames(gc), function(x) {
            dt <- as.data.table(t((t(fc) - gc[, x]) ^ 2 / gc[, x]))
            dt[which(is.na(dt), arr.ind = TRUE)] <- 0 # gc > 0
            dt[cTobject@counts <= 0] <- 0 # counts > 0
            ba <- vapply(gCobject@cl, function(y)
                dt[, Reduce('+', .SD), .SDcols = y],
                numeric(length = nrow(dt)))
            rowSums(ba[, deg > 1] * log10(csums[, deg > 1]),
                    na.rm = TRUE) / rowSums(A)
        }, numeric(length = nrow(cTobject@counts)))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ENC' Npvembre 2002.
###
#' @rdname codonUsage
#' @export
setGeneric(
    name = "ENCprime",
    def = function(cTobject,
                    subsets = list(),
                    self = TRUE,
                    ribosomal = FALSE,
                    id_or_name2 = "1",
                    alt.init = TRUE,
                    stop.rm = TRUE,
                    filtering = "none",
                    len.threshold = 80)
    {
        standardGeneric("ENCprime")
    }
)

#' @rdname codonUsage
setMethod(
    f = "ENCprime",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject,
                            subsets,
                            self,
                            ribosomal,
                            id_or_name2,
                            alt.init,
                            stop.rm,
                            filtering,
                            len.threshold)
    {
        if (filtering == "hard") {
            cTobject <- subset(cTobject, cTobject@len > len.threshold)
        } else if (filtering == "soft") {
            if (any(cTobject@len < len.threshold))
                warning("Some sequences have below-threshold length!")
        } else if (filtering == "none")
            NULL

        gCobject <- genCode(id_or_name2, alt.init, stop.rm)
        if (stop.rm) {
            cTobject@counts <- cTobject@counts[, gCobject@nostops]
            cTobject@len <-
                rowSums(cTobject@counts[, gCobject@nostops])
        }

        deg <- gCobject@deg
        csums <- .countsbyaa(cTobject, gCobject)
        fc <- .normFrequencies(cTobject, gCobject)
        gc <-
            expectedCU(cTobject, gCobject, subsets, self, ribosomal)
        vapply(colnames(gc), function(x) {
            # chi squared
            dt <- as.data.table(t((t(fc) - gc[, x]) ^ 2 / gc[, x]))
            dt[which(is.na(dt), arr.ind = TRUE)] <- 0
            chisum <- vapply(gCobject@cl, function(x)
                dt[, Reduce('+', .SD), .SDcols = x],
                numeric(length = nrow(dt)))
            chisq <- csums * chisum
            # homozygosity
            fa <- as.data.table(t((t(
                chisq + csums
            ) - deg) /
                (deg * t(csums - 1))))
            fa[csums < 5] <- NA
            .effNc(fa, deg)
        }, numeric(length = nrow(cTobject@counts)))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ENC, Wright 1990.
###
#' @rdname codonUsage
#' @export
setGeneric(
    name = "ENC",
    def = function(cTobject,
                    id_or_name2 = "1",
                    alt.init = TRUE,
                    stop.rm = TRUE,
                    filtering = "none",
                    len.threshold = 80)
    {
        standardGeneric("ENC")
    }
)

#' @rdname codonUsage
setMethod(
    f = "ENC",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject,
                            id_or_name2,
                            alt.init,
                            stop.rm,
                            filtering,
                            len.threshold)
    {
        if (filtering == "hard") {
            cTobject <- subset(cTobject, cTobject@len > len.threshold)
        } else if (filtering == "soft") {
            if (any(cTobject@len < len.threshold))
                warning("Some sequences have below-threshold length!")
        } else if (filtering == "none")
            NULL

        gCobject <- genCode(id_or_name2, alt.init, stop.rm)
        if (stop.rm) {
            cTobject@counts <- cTobject@counts[, gCobject@nostops]
            cTobject@len <-
                rowSums(cTobject@counts[, gCobject@nostops])
        }
        cl <- gCobject@cl
        counts <- as.data.table(cTobject@counts)
        csums <- .countsbyaa(cTobject, gCobject)
        freqs <- lapply(seq_along(cl), function(x)
            counts[, .SD / csums[, x], .SDcols = cl[[x]]])
        pi <- vapply(freqs, function(x)
            rowSums(x ^ 2, na.rm = TRUE),
            numeric(length = nrow(counts)))
        # remove AA with count <= 1
        csums[which(csums <= 1, arr.ind = TRUE)] <- NA
        fa <- as.data.table((csums * pi - 1) / (csums - 1))
        .effNc(fa, gCobject@deg)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SCUO, Wan et al. 2004.
###
#' @rdname codonUsage
#' @export
setGeneric(
    name = "SCUO",
    def = function(cTobject,
                    id_or_name2 = "1",
                    alt.init = TRUE,
                    stop.rm = FALSE,
                    filtering = "none",
                    len.threshold = 80)
    {
        standardGeneric("SCUO")
    }
)

#' @rdname codonUsage
setMethod(
    f = "SCUO",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject,
                            id_or_name2,
                            alt.init,
                            stop.rm,
                            filtering,
                            len.threshold)
    {
        if (filtering == "hard") {
            cTobject <- subset(cTobject, cTobject@len > len.threshold)
        } else if (filtering == "soft") {
            if (any(cTobject@len < len.threshold))
                warning("Some sequences have below-threshold length!")
        } else if (filtering == "none")
            NULL

        gCobject <- genCode(id_or_name2, alt.init, stop.rm)
        if (stop.rm) {
            cTobject@counts <- cTobject@counts[, gCobject@nostops]
            cTobject@len <-
                rowSums(cTobject@counts[, gCobject@nostops])
        }

        cl <- gCobject@cl
        counts <- as.data.table(cTobject@counts)
        csums <- .countsbyaa(cTobject, gCobject)
        freqs <- lapply(seq_along(cl), function(x)
            counts[, .SD / csums[, x], .SDcols = cl[[x]]])
        Ha <-
            vapply(freqs, function(x)
                rowSums(-x * log10(x), na.rm = TRUE),
                numeric(length = nrow(counts)))
        Hmax <- log10(vapply(cl, length, numeric(length = 1)))
        Oa <- t((Hmax - t(Ha)) / Hmax)
        Fa <- csums / rowSums(csums[, gCobject@deg > 1])
        rowSums(Oa * Fa, na.rm = TRUE)
    }
)
