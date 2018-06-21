#' @include codonTable-class.R
#' @include genCode-class.R
#' @include functions.R
#' @include codonUsage.R
#' @import data.table
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MELP, Supek and Vlahovicek 2005.
###
#' Calculate MILC-based Expression Level Predictor (MELP).
#'
#' Calculate MELP values for every sequence in the given
#' \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with MELP values for every specified subset
#'    (\code{subsets}, \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of MELP, see
#'    \href{https://bit.ly/2GkT7qe}{Supek & Vlahovicek (2005)}.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#'
#' # calculate MELP with respect to the CU
#' # of ribosomal genes among the example DNA sequences
#' melp <- MELP(cT, ribosomal = TRUE)
#' head(melp)
#'
#' # calculate MELP distance with respect to the average CU
#' # of the first 20 example DNA sequences
#' # (i.e. the first half of the example DNA set)
#' melp <- MELP(cT, subsets = list(half = c(rep(TRUE, 20), rep(FALSE, 20))))
#'
#' # alternatively, you can specify codonTable as a subset
#' halfcT <- codonTable(codonCounts(cT)[1:20,])
#' melp2 <- MELP(cT, subsets = list(half = halfcT))
#' all.equal(melp, melp2) # TRUE
#'
#' # filtering
#' MELP(cT, ribosomal = TRUE,
#'      filtering = "hard", len.threshold = 80) # MELP for 9 sequences
#'                                              # (note that, accidentally,
#'                                              # all are ribosomal)
#' sum(getlen(cT) > 80) # 9 sequences are longer than 80 codons
#' melp1 <- MELP(cT, ribosomal = TRUE, filtering = "none") # no filtering
#' melp2 <- MELP(cT, ribosomal = TRUE, filtering = "soft") # warning
#' all.equal(melp1, melp2) # TRUE
#'
#' # options for genetic code
#' melp <- MELP(cT, ribosomal = TRUE,
#'              stop.rm = TRUE) # don't use stop codons in calculation
#' melp <- MELP(cT, ribosomal = TRUE,
#'              alt.init = FALSE) # don't use alternative start codons
#' melp <- MELP(cT, ribosomal = TRUE,
#'              id_or_name2 = "2") # use different genetic code, for help
#'                                 # see `?Biostrings::GENETIC_CODE`
#'
#' @rdname MELP
#' @export
setGeneric(
    name = "MELP",
    def = function(cTobject,
                   subsets = list(),
                   ribosomal = FALSE,
                   id_or_name2 = "1",
                   alt.init = TRUE,
                   stop.rm = FALSE,
                   filtering = "none",
                   len.threshold = 80) {
        standardGeneric("MELP")
    }
)

#' @rdname MELP
#' @export
setMethod(
    f = "MELP",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject,
                          subsets,
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

        milcs <- MILC(cTobject,
                      subsets,
                      self = TRUE,
                      ribosomal,
                      id_or_name2,
                      alt.init,
                      stop.rm)
        vapply(as.list(colnames(milcs)[-1]), function(y)
            milcs[, "self"] / milcs[, y],
            numeric(length = nrow(milcs)))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### E, Karlin and Mrazek 2000
###
#' Calculate gene expression measure (E).
#'
#' Calculate E values for every sequence in the given
#' \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with E values for every specified subset (\code{subsets},
#'    \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of gene expression measure (E), see
#'    \href{https://bit.ly/2pGLno5}{Karlin and Mrazek (2000)}.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#'
#' # calculate E with respect to the CU
#' # of ribosomal genes among the example DNA sequences
#' em <- E(cT, ribosomal = TRUE)
#' head(em)
#'
#' # calculate E distance with respect to the average CU
#' # of the first 20 example DNA sequences
#' # (i.e. the first half of the example DNA set)
#' em <- E(cT, subsets = list(half = c(rep(TRUE, 20), rep(FALSE, 20))))
#'
#' # alternatively, you can specify codonTable as a subset
#' halfcT <- codonTable(codonCounts(cT)[1:20,])
#' em2 <- E(cT, subsets = list(half = halfcT))
#' all.equal(em, em2) # TRUE
#'
#' # filtering
#' E(cT, ribosomal = TRUE,
#'   filtering = "hard", len.threshold = 80) # E for 9 sequences (note that,
#'                                           # accidentally, all are ribosomal)
#' sum(getlen(cT) > 80) # 9 sequences are longer than 80 codons
#' em1 <- E(cT, ribosomal = TRUE, filtering = "none") # no filtering
#' em2 <- E(cT, ribosomal = TRUE, filtering = "soft") # warning
#' all.equal(em1, em2) # TRUE
#'
#' # options for genetic code
#' em <- E(cT, ribosomal = TRUE,
#'         stop.rm = TRUE) # don't use stop codons in calculation
#' em <- E(cT, ribosomal = TRUE,
#'         alt.init = FALSE) # don't use alternative start codons
#' em <- E(cT, ribosomal = TRUE,
#'         id_or_name2 = "2") # use different genetic code, for help
#'                            # see `?Biostrings::GENETIC_CODE`
#'
#' @rdname E
#' @export
setGeneric(
    name = "E",
    def = function(cTobject,
                   subsets = list(),
                   ribosomal = FALSE,
                   id_or_name2 = "1",
                   alt.init = TRUE,
                   stop.rm = FALSE,
                   filtering = "none",
                   len.threshold = 80) {
        standardGeneric("E")
    }
)

#' @rdname E
#' @export
setMethod(
    f = "E",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject,
                          subsets,
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

        b <- B(cTobject,
               subsets,
               self = TRUE,
               ribosomal,
               id_or_name2,
               alt.init,
               stop.rm)
        vapply(colnames(b)[-1], function(y)
            b[, "self"] / b[, y],
            numeric(length = nrow(b)))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### CAI, Sharp & Lee 1987
###
#' Calculate Codon Adaptation Index (CAI).
#'
#' Calculate CAI values for every sequence in the given
#' \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with CAI values for every specified subset
#'    (\code{subsets}, \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of Codon Adaptation Index (CAI), see
#'    \href{https://bit.ly/2HZuk8n}{Sharp and Li (1987)}.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#'
#' # calculate CAI with respect to the CU
#' # of ribosomal genes among the example DNA sequences
#' cai <- CAI(cT, ribosomal = TRUE)
#' head(cai)
#'
#' # calculate CAI distance with respect to the average CU
#' # of the first 20 example DNA sequences
#' # (i.e. the first half of the example DNA set)
#' cai <- CAI(cT, subsets = list(half = c(rep(TRUE, 20), rep(FALSE, 20))))
#'
#' # alternatively, you can specify codonTable as a subset
#' halfcT <- codonTable(codonCounts(cT)[1:20,])
#' cai2 <- CAI(cT, subsets = list(half = halfcT))
#' all.equal(cai, cai2) # TRUE
#'
#' # filtering
#' CAI(cT, ribosomal = TRUE,
#'   filtering = "hard", len.threshold = 80) # CAI for 9 sequences (note that,
#'                                           # accidentally, all are ribosomal)
#' sum(getlen(cT) > 80) # 9 sequences are longer than 80 codons
#' cai1 <- CAI(cT, ribosomal = TRUE, filtering = "none") # no filtering
#' cai2 <- CAI(cT, ribosomal = TRUE, filtering = "soft") # warning
#' all.equal(cai1, cai2) # TRUE
#'
#' # options for genetic code
#' cai <- CAI(cT, ribosomal = TRUE,
#'         stop.rm = TRUE) # don't use stop codons in calculation
#' cai <- CAI(cT, ribosomal = TRUE,
#'         alt.init = FALSE) # don't use alternative start codons
#' cai <- CAI(cT, ribosomal = TRUE,
#'         id_or_name2 = "2") # use different genetic code, for help
#'                            # see `?Biostrings::GENETIC_CODE`
#'
#' @rdname CAI
#' @export
setGeneric(
    name = "CAI",
    def = function(cTobject,
                   subsets = list(),
                   ribosomal = FALSE,
                   id_or_name2 = "1",
                   alt.init = TRUE,
                   stop.rm = FALSE,
                   filtering = "none",
                   len.threshold = 80) {
        standardGeneric("CAI")
    }
)

#' @rdname CAI
#' @export
setMethod(
    f = "CAI",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject,
                          subsets,
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
        deg <- gCobject@deg
        aa <- unlist(gCobject@ctab[, "AA"], use.names = FALSE)
        counts <- as.data.table(cTobject@counts)
        csums <- .countsbyaa(cTobject, gCobject)
        fc <- .normFrequencies(cTobject, gCobject)
        gc <-
            expectedCU(cTobject, gCobject, subsets, self = FALSE, ribosomal)
        gcmax <-
            apply(gc, 2, function(x)
                # freq of the most frequent codons by aa
                by(x, droplevels(aa), max))
        gcmax[gcmax == 0] <-
            0.5 # account for codons not used in ref. set
        x <- as.integer(droplevels(aa))
        gcm <-
            vapply(colnames(gc), function(y)
                # freq of the most freq. codons
                gc[, y] <- gcmax[x, y],
                numeric(length = nrow(gc)))
        nondeg <- match(which(deg == 1), as.integer(aa))
        counts[, (nondeg) := NA] # remove codons for non-degenerate aa
        vapply(colnames(gcm), function(y)
            exp(
                rowSums(counts * t(log(t(
                    fc
                ) / gcm[, y])), na.rm = TRUE) /
                    rowSums(counts, na.rm = TRUE)
            ),
            numeric(length = nrow(counts)))

    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Fop, Ikemura 1981
###
#' Calculate frequency of optimal codons (Fop).
#'
#' Calculate Fop values for every sequence in the given
#' \code{codonTable} object.
#'
#' @inheritParams MILC
#'
#' @return A matrix with Fop values for every specified subset
#'    (\code{subsets}, \code{self}, \code{ribosomal}) in columns.
#'
#'    For definition of frequency of optimal codons (Fop), see
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/6175758}{Ikemura (1981)}.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#'
#' # calculate Fop with respect to the CU
#' # of ribosomal genes among the example DNA sequences
#' fop <- Fop(cT, ribosomal = TRUE)
#' head(fop)
#'
#' # calculate Fop distance with respect to the average CU
#' # of the first 20 example DNA sequences
#' # (i.e. the first half of the example DNA set)
#' fop <- Fop(cT, subsets = list(half = c(rep(TRUE, 20), rep(FALSE, 20))))
#'
#' # alternatively, you can specify codonTable as a subset
#' halfcT <- codonTable(codonCounts(cT)[1:20,])
#' fop2 <- Fop(cT, subsets = list(half = halfcT))
#' all.equal(fop, fop2) # TRUE
#'
#' # filtering
#' Fop(cT, ribosomal = TRUE,
#'   filtering = "hard", len.threshold = 80) # Fop for 9 sequences
#' sum(getlen(cT) > 80) # 9 sequences are longer than 80 codons
#' fop1 <- Fop(cT, ribosomal = TRUE, filtering = "none") # no filtering
#' fop2 <- Fop(cT, ribosomal = TRUE, filtering = "soft") # warning
#' all.equal(fop1, fop2) # TRUE
#'
#' # options for genetic code
#' fop <- Fop(cT, ribosomal = TRUE,
#'            stop.rm = TRUE) # don't use stop codons in calculation
#' fop <- Fop(cT, ribosomal = TRUE,
#'            alt.init = FALSE) # don't use alternative start codons
#' fop <- Fop(cT, ribosomal = TRUE,
#'            id_or_name2 = "2") # use different genetic code, for help
#'                               # see `?Biostrings::GENETIC_CODE`
#'
#' @rdname Fop
#' @export
setGeneric(
    name = "Fop",
    def = function(cTobject,
                   subsets = list(),
                   ribosomal = FALSE,
                   id_or_name2 = "1",
                   alt.init = TRUE,
                   stop.rm = FALSE,
                   filtering = "none",
                   len.threshold = 80) {
        standardGeneric("Fop")
    }
)

#' @rdname Fop
#' @export
setMethod(
    f = "Fop",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject,
                          subsets,
                          ribosomal,
                          id_or_name2,
                          alt.init,
                          stop.rm,
                          filtering,
                          len.threshold)  {
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
        aa <- unlist(gCobject@ctab[, "AA"], use.names = FALSE)
        counts <- as.data.table(cTobject@counts)
        csums <- .countsbyaa(cTobject, gCobject)
        fc <- .normFrequencies(cTobject, gCobject)
        gc <-
            expectedCU(cTobject, gCobject, subsets, self = FALSE, ribosomal)
        gcmax <-
            apply(gc, 2, function(x)
                # freq of the most frequent codons by aa
                by(x, droplevels(aa), max))
        gcmax[gcmax == 0] <-
            0.5 # account for codons not used in ref. set
        x <- as.integer(droplevels(aa))
        gcm <-
            vapply(colnames(gc), function(y)
                # freq of the most freq. codons
                gc[, y] <- gcmax[x, y], numeric(length = nrow(gc)))
        nondeg <- match(which(deg == 1), as.integer(aa))
        counts[, (nondeg) := NA] # remove codons for non-degenerate aa
        vapply(colnames(gcm), function(y) {
            ra <- t(fc) / gcm[, y] # relative adaptivenes of codons
            cnt <- as.matrix(counts)
            cnt[which(t(ra) < 0.9, arr.ind = TRUE)] <-
                NA # as impl. in INCA
            rowSums(cnt, na.rm = TRUE) / rowSums(counts, na.rm = TRUE)
        }, numeric(length = nrow(counts)))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GCB, Merkl 2003
###
#' Calculate gene codon bias (GCB).
#'
#' Calculate GCB values for every sequence in the given
#' \code{codonTable} object.
#'
#' @inheritParams MILC
#' @param seed A logical vector, of the length equal to
#'    \code{getlen(cTobject)}, or a character vector (of any length)
#'    containing KEGG/eggNOG annotations, or a codonTable object
#'    (of any length). Indicates a set of genes, or their CU, to be used
#'    as a target in the first iteration of the algorithm.
#' @param ribosomal Logical, if \code{TRUE}, ribosomal genes are used
#'    as a seed. Default is \code{FALSE}, in which case \code{seed}
#'    has to be specified.
#' @param perc percent of top ranking genes to be used as a target set
#'    for the next iteration. Default is 0.05.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#'
#' # calculate GCB with CU of ribosomal genes among the example DNA sequences
#' # used as a target (seed) in the first iteration of the algorithm
#' gcb <- GCB(cT, ribosomal = TRUE)
#' head(gcb)
#'
#' # calculate GCB distance with the first 20 example DNA sequences
#' # (i.e. the first half of the example DNA set) as a seed
#' gcb <- GCB(cT, seed = c(rep(TRUE, 20), rep(FALSE, 20)))
#'
#' # alternatively, you can specify codonTable as a seed
#' halfcT <- codonTable(codonCounts(cT)[1:20,])
#' gcb2 <- GCB(cT, seed = halfcT)
#' all.equal(gcb, gcb2) # TRUE
#'
#' # options for genetic code
#' gcb <- GCB(cT, ribosomal = TRUE,
#'            stop.rm = TRUE) # don't use stop codons in calculation
#' gcb <- GCB(cT, ribosomal = TRUE,
#'            alt.init = FALSE) # don't use alternative start codons
#' gcb <- GCB(cT, ribosomal = TRUE,
#'            id_or_name2 = "2") # use different genetic code, for help
#'                               # see `?Biostrings::GENETIC_CODE`
#'
#' @return A numeric vector with GCB values.
#'
#'    For definition of gene codon bias (GCB), see
#'    \href{http://www.ncbi.nlm.nih.gov/pubmed/14708578}{Merkl (2003)}.
#'

#' @rdname GCB
#' @export
setGeneric(
    name = "GCB",
    def = function(cTobject,
                   seed = logical(),
                   ribosomal = FALSE,
                   perc = 0.05,
                   id_or_name2 = "1",
                   alt.init = TRUE,
                   stop.rm = FALSE,
                   filtering = "none",
                   len.threshold = 80) {
        standardGeneric("GCB")
    }
)

#' @rdname GCB
#' @export
setMethod(
    f = "GCB",
    signature = c(cTobject = "codonTable"),
    definition = function(cTobject,
                          seed,
                          ribosomal,
                          perc,
                          id_or_name2,
                          alt.init,
                          stop.rm,
                          filtering,
                          len.threshold)  {
        if (filtering == "hard") {
            cTobject <- subset(cTobject, cTobject@len > len.threshold)
        } else if (filtering == "soft") {
            if (any(cTobject@len < len.threshold))
                warning("Some sequences have below-threshold length!")
        } else if (filtering == "none")
            NULL

        if (length(seed) == 0 & ribosomal == FALSE)
            stop("Seed is not specified!")

        gCobject <- genCode(id_or_name2, alt.init, stop.rm)
        if (stop.rm) {
            cTobject@counts <- cTobject@counts[, gCobject@nostops]
            cTobject@len <-
                rowSums(cTobject@counts[, gCobject@nostops])
        }

        len <- cTobject@len
        counts <- as.data.table(cTobject@counts)
        fc <- .normFrequencies(cTobject, gCobject)
        if (length(seed) == 0)
            seed <- list()
        else
            seed <- list(seed = seed)
        gc <-
            expectedCU(cTobject, gCobject, seed, self = FALSE, ribosomal)
        gcb_prev <- numeric(nrow(counts))
        iter <- 0
        repeat {
            cb <- log(gc / colMeans(fc, na.rm = TRUE))
            cb[gc == 0] <- -5
            gcb <-
                rowSums(t(t(counts) * as.numeric(cb)), na.rm = TRUE) / len
            diff <- all(gcb == gcb_prev)
            if (diff | iter > 6)
                break
            else {
                iter <- iter + 1
                gcb_prev <- gcb
                top <-
                    order(gcb, decreasing = TRUE)[seq_len(perc * length(gcb))]
                gc <- expectedCU(cTobject,
                                 gCobject,
                                 list(s = codonTable(counts[top, ])),
                                 self = FALSE,
                                 ribosomal)
            }
        }
        return(gcb)
    }
)
