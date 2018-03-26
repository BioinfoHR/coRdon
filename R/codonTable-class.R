#' @include genCode-class.R
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom Biostrings width
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
#' @slot COG  A character vector of COG annotations for sequences, length
#'    equal to \code{nrow(counts)}. If no annotation is available, this
#'    will be an empty vector.
#'
#' @examples
#' # create codonTable with codon counts for sequences in DNAStringSet
#' require(Biostrings)
#' dna <- DNAStringSet(c("ACGAAGTGTACTGTAATTTGCACAGTACTTAAATGT",
#'                       "ACGTCCGTACTGATCGATTCCGTGATT"))
#' cT <- codonTable(dna)
#' codonCounts(cT)
#' getlen(cT)
#' getKO(cT)
#' cT <- setKO(cT, c("K00001", "K00002"))
#' getKO(cT)
#'
#' # convert matrix containing codon counts to codonTable
#' mat <- matrix(sample(1:10, 122, replace = TRUE), nrow = 2)
#' codonTable(mat) # produces informative warning
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
            msg <- paste("Codon counts have to be integers. \n")
            errors <- c(errors, msg)
        }
        if (all(rowSums(object@counts) != object@len)) {
            msg <-
                "(Some of) summed codon counts differ from sequence length.\n"
            errors <- c(errors, msg)
        }
        if (KOlen != 0 & KOlen != ns) {
            msg <- cat("Number of KO annotations,", KOlen,
                       "differ from the number of sequences,", ns, "\n")
            errors <- c(errors, msg)
        }
        if (COGlen != 0 & COGlen != ns) {
            msg <- cat("Number of COG annotations,", COGlen,
                       " differ from the number of sequences,", ns, "\n")
            errors <- c(errors, msg)
        }
        if (length(errors) == 0) TRUE else stop(errors)
    }
)

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
#' @param x An object of \code{DNAStringSet}, \code{matrix} or
#'   \code{data.frame} class.
#'
#' @return A \code{codonTable}.
#'
#' @export
setMethod(
    f = "codonTable",
    signature = "DNAStringSet",
    definition = function(x) {
        bad <- which(width(x) %% 3 != 0)
        if (length(bad) != 0)
            warning(paste0(
                "\nLength of sequence(s) at the following postion(s) ",
                "is not divisible by 3: \n",
                bad,
                ".\nDiscarding surplus nucleotides.\n"))
        ctb <- .codonTable(x)
        ctb <- ctb[,order(colnames(ctb))] # sort codons alphabetically
        if (class(ctb) == "integer") # in case there is only one sequence
            ctb <- rbind(ctb)
        new(
            "codonTable",
            ID = if (!is.null(names(x))) names(x) else character(),
            counts = ctb,
            len = unname(rowSums(ctb)),
            KO = regmatches(names(x),
                            regexpr("K\\d{5}", names(x))),
            COG = regmatches(names(x),
                             regexpr("([KCN]|TW)OG\\d{5}", names(x)))
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
            if (is.null(colnames(x)) |
                !all(sort(colnames(x))==genCode()@codons)) {
                colnames(x) <- genCode()@codons
                message(
                    "Assigned alphabetically sorted codons as column names."
                    )
            }
        } else if (ncol(x) == 61) {
            if (is.null(colnames(x))) {
                colnames(x) <- genCode()@nostops
                message(
                    "Assigned alphabetically sorted codons as column names."
                    )
            } else if (!all.equal(sort(colnames(x)), genCode()@nostops)) {
                colnames(x) <- genCode()@nostops
                message(
                    "Assigned alphabetically sorted codons as column names."
                    )
            }
        } else stop("Matrix should have 64 or 61 columns!")

        new(
            "codonTable",
            ID = if (!is.null(rownames(x))) rownames(x) else character(),
            counts = rbind(x[,sort(colnames(x))]), # sort codons alphabet.
            len = rowSums(as.matrix(x), na.rm = TRUE),
            KO = regmatches(rownames(x),
                            regexpr("K\\d{5}", rownames(x))),
            COG = regmatches(rownames(x),
                             regexpr("([KCN]|TW)OG\\d{5}", rownames(x)))
        )
    }
)

#' @rdname codonTable-class
#' @export
setMethod(
    f = "codonTable",
    signature = "data.frame",
    definition = function(x){
        x <- as.matrix(x)
        codonTable(x)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### codonTable accesor methods
###

#' Display the object of \code{codonTable} class.
#'
#' @docType methods
#' @name show-codonTable
#' @rdname show-codonTable
#' @aliases show-codonTable show,codonTable-method
#'
#' @param object A \code{codonTable} object.
#'
#' @return \code{show} returns an invisible \code{NULL}.
#'
#' @export
setMethod(
    f = "show",
    signature = "codonTable",
    definition = function(object){
        ns <- nrow(object@counts)
        IDs <- capture.output(str(object@ID))
        lens <- capture.output(str(object@len))
        KOs <- capture.output(str(object@KO))
        COGs <- capture.output(str(object@COG))
        cat("codonTable instance with codon counts from", ns, "sequences.\n")
        if (length(object@ID) != 0 & any(!is.na(object@ID)))
            cat("sequence IDs:\n", IDs, "\n")
        cat("sequence lengths:\n", lens, "\n")
        if (length(object@KO) != 0 & any(!is.na(object@KO)))
            cat("KO annotations:\n", KOs, "\n")
        if (length(object@COG) != 0 & any(!is.na(object@COG)))
            cat("COG annotations:\n", COGs, "\n")
    }
)

#' Length of \code{codonTable} object.
#'
#' Length of \code{codonTable} object is the number of sequences
#' for which there are codon counts contained in the object.
#'
#' @docType methods
#' @name length-codonTable
#' @rdname length-codonTable
#' @aliases length-codonTable length,codonTable-method
#'
#' @param x A \code{codonTable} object.
#'
#' @return Numeric, the length of \code{x}.
#'
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
    name = "codonCounts",
    def = function(object){
        standardGeneric("codonCounts")
    }
)

#' @describeIn codonTable Get codon counts from \code{codonTable} object.
#'
#' @param object A \code{codonTable} object.
#'
#' @export
setMethod(
    f = "codonCounts",
    signature = "codonTable",
    definition = function(object){
        return(object@counts)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "getID",
    def = function(object){
        standardGeneric("getID")
    }
)

#' @describeIn codonTable Get IDs for \code{codonTable} class.
#'
#' @inheritParams codonCounts
#'
#' @export
setMethod(
    f = "getID",
    signature = "codonTable",
    definition = function(object){
        return(object@ID)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "getlen",
    def = function(object){
        standardGeneric("getlen")
    }
)
#' @describeIn codonTable
#' Get lengths of sequences in \code{codonTable} object.
#'
#' @inheritParams codonCounts
#'
#' @export
setMethod(
    f = "getlen",
    signature = "codonTable",
    definition = function(object){
        return(object@len)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "getKO",
    def = function(object){
        standardGeneric("getKO")
    }
)
#' @describeIn codonTable
#' Get KO annotations of sequences in \code{codonTable} object.
#'
#' @inheritParams codonCounts
#'
#' @export
setMethod(
    f = "getKO",
    signature = "codonTable",
    definition = function(object){
        return(object@KO)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "setKO",
    def = function(object, ann){
        standardGeneric("setKO")
    }
)
#' @describeIn codonTable Set KO annotations for \code{codonTable} object.
#'
#' @inheritParams codonCounts
#' @param ann A character vector of sequence annotations, must be of length
#'    equal to \code{length(object)}.
#'
#' @export
setMethod(
    f = "setKO",
    signature = "codonTable",
    definition = function(object, ann){
        object@KO <- ann
        validObject(object)
        return(object)
    }
)

#' @rdname codonTable-class
#' @export
setGeneric(
    name = "getCOG",
    def = function(object){
        standardGeneric("getCOG")
    }
)

#' @describeIn codonTable
#' Get COG annotations of sequences in \code{codonTable} object.
#'
#' @inheritParams codonCounts
#'
#' @export
setMethod(
    f = "getCOG",
    signature = "codonTable",
    definition = function(object){
        return(object@COG)
    }
)


#' @rdname codonTable-class
#' @export
setGeneric(
    name = "setCOG",
    def = function(object, ann){
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
    definition = function(object, ann){
        object@COG <- ann
        validObject(object)
        return(object)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### codonTable subset methods
###
setGeneric(
    name = "mySubset",
    def = function(x, subset){
        stop("this is only a generic function: it should never be called!")
    }
)

setMethod(
    f = "mySubset",
    signature = c(x = "codonTable", subset = "logical"),
    definition = function(x, subset){
        if (length(x@ID[subset]) == 0) stop("Empty codonTable object!")
        nms <- getID(x)[subset]
        new("codonTable",
            ID = if(all(is.na(nms))) character() else nms,
            counts = rbind(x@counts[subset, ]),
            len = x@len[subset],
            KO = x@KO[subset],
            COG = x@COG[subset])
    }
)

setMethod(
    f = "mySubset",
    signature = c(x = "codonTable", subset = "character"),
    definition = function(x, subset){
        KOs <- x@KO %in% subset
        COGs <- x@COG %in% subset
        if (any(KOs)) {
            s <- KOs
        } else if (any(COGs)) {
            s <- COGs
        } else stop("No sequence has given annotation!")
        nms <- x@ID[s]
        new("codonTable",
            ID = if(all(is.na(nms))) character() else nms,
            counts = rbind(x@counts[s, ]),
            len = x@len[s],
            KO = x@KO[s],
            COG = x@COG[s])
    }
)

setMethod(
    f = "mySubset",
    signature = c(x = "codonTable", subset = "ANY"),
    definition = function(x, subset){
        stop("Object to be subset should be of class codonTable!")
    }
)

#' Subset \code{codonTable} object.
#'
#' @param x A \code{codonTable} object to be subset.
#' @param subset A logical or character vector indicating which elements of
#'    \code{x} to keep. If logical, \code{subset} should be of length
#'    \code{length(x)}. If character, \code{subset} should contain
#'     at least some of the elements of either \code{getKO(x)} or
#'     \code{getCOG(x)}.
#' @inheritParams base::subset
#'
#' @return subsets of \code{codonTable} object, keeping in each slot only
#' those elements that meet the criteria in \code{subset}.
#'
#' @examples
#' # create codonTable
#' mat <- matrix(sample(1:10, 610, replace = TRUE), nrow = 10)
#' cT <- codonTable(mat) # produces informative warning
#' cT
#' subset(cT, c(rep(c(TRUE,FALSE), 5))) # subset odd sequences
#'
#' cT <- setKO(cT, rep(c("K00001", "K00002"), 5))
#' subset(cT, "K00001")
#'
#' cT <- setCOG(cT, rep(c("COG0001", "COG0002"), 5))
#' subset(cT, "COG0001")
#'
#' @rdname subset
#' @export
subset.codonTable <- function(x, subset, ...){
    mySubset(x,subset)
}
