#' An S4 class \code{genCode}
#'
#' Object of \code{genCode} class describes the variant of genetic code
#' to be used in CU calculations.
#'
#' @slot ctab A \code{data.table} with two colums: \code{codon} and \code{AA},
#'    amino acid.
#' @slot codons A character vector of codons.
#' @slot stops A character vector of stop codons. Note that, if \code{stop.rm} is
#'    \code{TRUE}, this will be an empty vector.
#' @slot nostops A character vector of no-stop codons. If \code{stop.rm} is
#'    \code{TRUE}, this will be equal to the \code{codons} slot.
#' @slot cl A list, each element of which is a vector of integers indicating
#'    the positions of synonymous codons for that amino acid, when codons are
#'    alphabetically ordered.
#' @slot deg A numeric vector of degeneracies for alphabetically ordered amino
#'    acids.
#'
setClass(
    "genCode",
    slots = c(
        ctab = "data.table",
        codons = "character",
        stops = "character",
        nostops = "character",
        cl = "list",
        deg = "numeric"
    )
)

#' @importFrom Biostrings getGeneticCode
#' @import data.table
.genCode <- function(id_or_name2, alt.init, stop.rm)
{
    ctab <- as.data.table(getGeneticCode(id_or_name2, as.data.frame = T),
                          keep.rownames = "codon")
    ctab[, AA := as.factor(ctab$AA)]
    if (alt.init) ctab[Start == "M", AA:="M"]
    ctab[, Start := NULL]
    if (stop.rm) ctab <- ctab[AA != "*", ]
    return(ctab[order(ctab$codon),])
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### genCode constructor (not exported)
###

#' @rdname genCode-class
#' @importFrom Biostrings getGeneticCode
#' @import data.table
setGeneric(
    name = "genCode",
    def = function(id_or_name2 = "1", alt.init = TRUE, stop.rm = FALSE){
        standardGeneric("genCode")
    }
)

#' @describeIn genCode Creates new instances of \code{genCode} class.
#'
#' @inheritParams Biostrings::getGeneticCode
#' @param alt.init logical, whether to use alternative initiation codons.
#'   Default is \code{TRUE}.
#' @param stop.rm logical, whether to remove stop codons. Default is \code{FALSE}.
#'
setMethod(
    f = "genCode",
    definition = function(id_or_name2, alt.init, stop.rm){
        ctab <- .genCode(id_or_name2, alt.init, stop.rm)
        cl <- sapply(levels(droplevels(ctab$AA)),
                    function(x) which(ctab$AA==x), USE.NAMES = T)
        new("genCode",
            ctab = ctab,
            codons = ctab[, codon],
            stops = ctab[AA == "*", codon],
            nostops = ctab[AA != "*", codon],
            cl = cl,
            deg = sapply(cl, length))
    }
)
