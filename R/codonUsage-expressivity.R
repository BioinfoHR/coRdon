#' @include codonTable-class.R
#' @include genCode-class.R
#' @include functions.R
#' @include codonUsage.R
#' @import data.table
NULL

#' Calculate CU expressivity measures.
#'
#' Calculate values of the CU expressivity measure 
#' for every sequence in the given \code{codonTable} object.
#' The following methods are implemented:  
#'  \code{MELP}, CU expressivity measure based on
#'  Measure Independent of Length and Composition 
#'  \href{https://bit.ly/2GkT7qe}{Supek & Vlahovicek (2005)},
#'  \code{E}, gene expression measure (E) 
#'  \href{https://bit.ly/2pGLno5}{Karlin and Mrazek (2000)}, 
#'  \code{CAI}, Codon Adaptation Index (CAI) 
#'  \href{https://bit.ly/2HZuk8n}{Sharp and Li (1987)}, 
#'  \code{Fop}, frequency of optimal codons (Fop) 
#'  \href{https://www.ncbi.nlm.nih.gov/pubmed/6175758}{Ikemura (1981)}, 
#'  \code{GCB}, gene codon bias (GCB) 
#'  \href{http://www.ncbi.nlm.nih.gov/pubmed/14708578}{Merkl (2003)}.
#'  
#' @inheritParams codonUsage
#' @param seed A logical vector, of the length equal to
#'    \code{getlen(cTobject)}, or a character vector (of any length)
#'    containing KEGG/eggNOG annotations, or a codonTable object
#'    (of any length). Used only in GCB calculation. 
#'    Indicates a set of genes, or their CU, to be used
#'    as a target in the first iteration of the algorithm.
#' @param perc percent of top ranking genes to be used as a target set
#'    for the next iteration of the algorithm that calculates GCB. 
#'    Default is 0.05.
#'
#' @return A matrix (for GCB a numeric vector) with CU expressivity values 
#'    for every specified subset (\code{subsets}, \code{self}, 
#'    \code{ribosomal}) in columns.
#'
#' @examples
#' # load example DNA sequences
#' exampledir <- system.file("extdata", package = "coRdon")
#' cT <- codonTable(readSet(exampledir))
#' 
#' # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' # In the examples below, MELP values are calculated for all sequences; 
#' # any other CU expressivity measure can be caluclated in the same way,
#' # the only exception being GCB which takes `seed` instead of `subset` 
#' # parameter. (The exemples for GCB calculation are further below).
#' # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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
#' # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' # GCB calculationd
#' # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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
#' @name codonUsage-expressivity
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MELP, Supek and Vlahovicek 2005.
###
#' @rdname codonUsage-expressivity
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
                    len.threshold = 80)
    {
            standardGeneric("MELP")
    }
)

#' @rdname codonUsage-expressivity
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
                            len.threshold)
    {
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
        vapply(colnames(milcs)[-1], function(y)
            milcs[, "self"] / milcs[, y],
            numeric(length = nrow(milcs)))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### E, Karlin and Mrazek 2000
###
#' @rdname codonUsage-expressivity
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
                    len.threshold = 80)
    {
        standardGeneric("E")
    }
)

#' @rdname codonUsage-expressivity
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
                            len.threshold)
    {
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
#' @rdname codonUsage-expressivity
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
                    len.threshold = 80)
    {
        standardGeneric("CAI")
    }
)

#' @rdname codonUsage-expressivity
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
#' @rdname codonUsage-expressivity
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
                    len.threshold = 80)
    {
        standardGeneric("Fop")
    }
)

#' @rdname codonUsage-expressivity
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
#' @rdname codonUsage-expressivity
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
                    len.threshold = 80)
    {
        standardGeneric("GCB")
    }
)

#' @rdname codonUsage-expressivity
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
                            len.threshold)
    {
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
                gc <- expectedCU(
                        cTobject,
                        gCobject,
                        list(s = codonTable(counts[top, ])),
                        self = FALSE,
                        ribosomal)
            }
        }
        return(gcb)
    }
)
