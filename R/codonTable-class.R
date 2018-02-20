#' @include genCode-class.R
NULL

#' An S4 class \code{codonTable}
#'
#' Contains codon counts and optional annotation for a set DNA sequences.
#'
#' @slot ID A character vector of sequence identifiers.
#' @slot counts A matrix containing codon counts. Columns are codons, rows are
#'    sequences.
#' @slot len A numeric vector,length equal to \code{nrow(counts)}, containing
#'    lengths of sequnces.
#' @slot KO A character vector of KEGG annotations for sequences, length equal
#'    to \code{nrow(counts)}. If no annotation is available, this will be
#'    an empty vector.
#' @slot COG  A character vector of COG annotations for sequences, length equal
#'    to \code{nrow(counts)}. If no annotation is available, this will be
#'    an empty vector.
#'
setClass(
    "codonTable",
    slots = c(
        ID = "character",
        counts = "matrix",
        len = "numeric",
        KO = "character",
        COG = "character"
    )
)

setValidity(
    "codonTable",
    function(object){
        ns <- nrow(object@counts)
        KOlen <- length(object@KO)
        COGlen <- length(object@COG)
        errors <- character()
        if (!is.integer(object@counts)) {
            msg <- paste("Codon counts have to be integers!")
            errors <- c(errors, msg)
        }
        if (all(rowSums(object@counts) != object@len)) {
            msg <- "(Some of) summed codon counts differ from sequence length!"
            errors <- c(errors, msg)
        }
        if (KOlen != 0 & KOlen != ns) {
            msg <- cat("Number of KO annotations,", KOlen,
                       "differ from the number of sequences,", ns, "!")
            errors <- c(errors, msg)
        }
        if (COGlen != 0 & COGlen != ns) {
            msg <- cat("Number of COG annotations,", COGlen,
                       " differ from the number of sequences,", ns, "!")
            errors <- c(errors, msg)
        }
        if (length(errors) == 0) TRUE else errors
    }
)

#' @importFrom Biostrings oligonucleotideFrequency
.codonTable <- function(x)
{
    oligonucleotideFrequency(x, width = 3, step = 3)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### codonTable constructor
###

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "codonTable",
    def = function(x){
        standardGeneric("codonTable")
    }
)

#' @describeIn codonTable Create new objects of class \code{codonTable}.
#'
#' @param x An object of \code{DNAStringSet}, \code{matrix} or \code{data.frame} class.
#'
#' @export
setMethod(
    f = "codonTable",
    signature = "DNAStringSet",
    definition = function(x) {
        ctb <- .codonTable(x)
        new(
            "codonTable",
            ID = names(x),
            counts = ctb[,order(colnames(ctb))], # sort codons alphabetically
            len = rowSums(ctb),
            KO = regmatches(names(x), regexpr("K\\d{5}", names(x))),
            COG = regmatches(names(x), regexpr("([KCN]|TW)OG\\d{5}", names(x)))
        )
    }
)

#' @rdname codonTable-class
#' @export
setMethod(
    f = "codonTable",
    signature = "matrix",
    definition = function(x) {
        if (ncol(x) == 64) {
            if (is.null(colnames(x))) {
                colnames(x) <- genCode()@codons
                message("Assigned alphabetically sorted codons as column names.")
            }
        } else if (ncol(x) == 61) {
            if (is.null(colnames(x))) {
                colnames(x) <- genCode()@nostops
                message("Assigned alphabetically sorted codons as column names.")
            }
        } else stop("Matrix should have 64 or 61 columns!")

        new(
            "codonTable",
            ID = if (!is.null(rownames(x))) rownames(x) else character(),
            counts = rbind(x[,sort(colnames(x))]), # sort codons alphabetically
            len = rowSums(as.matrix(x), na.rm = TRUE),
            KO = regmatches(rownames(x), regexpr("K\\d{5}", rownames(x))),
            COG = regmatches(rownames(x), regexpr("([KCN]|TW)OG\\d{5}", rownames(x)))
        )
    }
)

#' @rdname codonTable-class
#' @export
setMethod(
    f = "codonTable",
    signature = "data.frame",
    definition = function(x)
        new(
            "codonTable",
            ID = if (!is.null(rownames(x))) rownames(x) else character(),
            counts = as.matrix(x)[,sort(colnames(x))], # sort codons alphabetically
            len = rowSums(x, na.rm = TRUE),
            KO = regmatches(rownames(x), regexpr("K\\d{5}", rownames(x))),
            COG = regmatches(rownames(x), regexpr("([KCN]|TW)OG\\d{5}", rownames(x)))
        )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### codonTable accesor methods
###

#' Display the object of \code{codonTable} class.
#'
#' @inheritParams methods::show
#'
#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = "codonTable",
    definition = function(object){
        ns <- nrow(object@counts)
        cat("codonTable instance with codon counts from", ns, "sequences")
    }
)

#' Length of \code{codonTable} object
#'
#' Get the length of \code{codonTable} object, i.e. the number of sequences
#' for which there are codon counts contained in the object.
#'
#' @param x A \code{codonTable} object.
#'
#' @rdname length
#' @export
setMethod(
    f = "length",
    signature = "codonTable",
    definition = function(x){
        nrow(codonCounts(x))
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "getID",
    def = function(cTobject){
        standardGeneric("getID")
    }
)

#' @describeIn codonTable Get IDs for \code{codonTable} class.
#'
#' @param cTobject A \code{codonTable} object.
#'
#' @export
setMethod(
    f = "getID",
    signature = "codonTable",
    definition = function(cTobject){
        return(cTobject@ID)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "codonCounts",
    def = function(cTobject){
        standardGeneric("codonCounts")
    }
)

#' @describeIn codonTable Get codon counts from \code{codonTable} object.
#'
#' @inheritParams getID
#'
#' @export
setMethod(
    f = "codonCounts",
    signature = "codonTable",
    definition = function(cTobject){
        return(cTobject@counts)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "getlen",
    def = function(cTobject){
        standardGeneric("getlen")
    }
)
#' @describeIn codonTable Get KO annotations for \code{codonTable} object.
#'
#' @inheritParams getID
#'
#' @export
setMethod(
    f = "getlen",
    signature = "codonTable",
    definition = function(cTobject){
        return(cTobject@len)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "getKO",
    def = function(cTobject){
        standardGeneric("getKO")
    }
)
#' @describeIn codonTable Get KO annotations for \code{codonTable} object.
#'
#' @inheritParams getID
#'
#' @export
setMethod(
    f = "getKO",
    signature = "codonTable",
    definition = function(cTobject){
        return(cTobject@KO)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "setKO",
    def = function(cTobject, ann){
        standardGeneric("setKO")
    }
)
#' @describeIn codonTable Set KO annotations for \code{codonTable} object.
#'
#' @inheritParams getID
#' @param ann A character vector of sequence annotations, must be of length
#'    equal to \code{length(cTobject)}.
#'
#' @export
setMethod(
    f = "setKO",
    signature = "codonTable",
    definition = function(cTobject, ann){
        cTobject@KO <- ann
        ok <- validObject(cTobject)
        if (ok) return(cTobject)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "getCOG",
    def = function(cTobject){
        standardGeneric("getCOG")
    }
)

#' @describeIn codonTable Get COG annotations for \code{codonTable} object.
#'
#' @inheritParams getID
#'
#' @export
setMethod(
    f = "getCOG",
    signature = "codonTable",
    definition = function(cTobject){
        return(cTobject@COG)
    }
)


#' @rdname codonTable-class
#' @export
setGeneric(
    name = "setCOG",
    def = function(cTobject, a){
        standardGeneric("setCOG")
    }
)

#' @describeIn codonTable Set COG annotations for \code{codonTable} object.
#'
#' @inheritParams setKO
#'
#' @export
setMethod(
    f = "setCOG",
    signature = "codonTable",
    definition = function(cTobject, a){
        cTobject@COG <- a
        ok <- validObject(cTobject)
        if (ok) return(cTobject)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### codonTable subset methods
###

#' Subsetting objects of \code{codonTable} class.
#'
#' Return subsets of \code{codonTable} object, keeping in each slot only those
#' elements that meet the criteria in \code{subset}.
#'
#' @rdname subset
#' @export
setGeneric(
    name = "subset",
    def = function(x, subset){
        standardGeneric("subset")
    }
)
#' @describeIn subset Subset \code{codonTable} object.
#'
#' @param x A \code{codonTable} object to be subset.
#' @param subset A logical or character vector indicating which elements of
#'    \code{cTobject} to keep. If logical, \code{subset} should be of length
#'    \code{nrow(counts(cTobject))}. If character, \code{subset} should contain
#'     at least some of the elements of either \code{getKO(codonTable)} or
#'     \code{getCOG(codonTable)}.
#'
#' @export
setMethod(
    f = "subset",
    signature = c("codonTable", "logical"),
    definition = function(x, subset){
        if (length(x@ID[subset]) == 0) stop("Empty codnoTable object!")
        new("codonTable",
            ID = getID(x)[subset],
            counts = rbind(x@counts[subset, ]),
            len = x@len[subset],
            KO = x@KO[subset],
            COG = x@COG[subset])
    }
)

#' @rdname subset
#' @export
setMethod(
    f = "subset",
    signature = c("codonTable", "character"),
    definition = function(x, subset){
        KOs <- x@KO %in% subset
        COGs <- x@COG %in% subset
        if (any(KOs)) {
            s <- KOs
        } else if (any(COGs)) {
            s <- COGs
        } else stop("No sequence has given annotation!")

        new("codonTable",
            ID = x@ID[s, ],
            counts = x@counts[s, ],
            len = x@len[s, ],
            KO = x@KO[s, ],
            COG = x@COG[s, ])
    }
)
