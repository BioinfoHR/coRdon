#' @include genCode-class.R
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom Biostrings width
#' @importFrom stringr str_extract
NULL

#' An S4 class \code{codonTable}
#'
#' Contains codon counts and optional annotation
#' for a set DNA sequences.
#'
#' @slot ID A character vector of sequence identifiers.
#' @slot counts A matrix containing codon counts.
#'    Columns are codons, rows are sequences.
#' @slot len A numeric vector,length equal to \code{nrow(counts)},
#'    containing lengths of sequnces.
#' @slot KO A character vector of KEGG annotations for sequences,
#'    length equal to \code{nrow(counts)}. If no annotation
#'    is available, this will be an empty vector.
#' @slot COG  A character vector of COG annotations for sequences,
#'    length equal to \code{nrow(counts)}. If no annotation
#'    is available, this will be an empty vector.
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
    warnings <- character()
    if (KOlen != 0 & KOlen != ns) {
      msg <- cat("Number of KO annotations,", KOlen,
                 "differ from the number of sequences,", ns, "\n")
      warnings <- c(warnings, msg)
    }
    if (COGlen != 0 & COGlen != ns) {
      msg <- cat("Number of COG annotations,", COGlen,
                 " differ from the number of sequences,", ns, "\n")
      warnings <- c(warnings, msg)
    }
    if (length(warnings) != 0) warning(warnings)
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
    empty <- which(width(x) == 0)
    if (length(empty) != 0) {
      warning(paste0(
        "\nLength of sequence(s) at the following postion(s) ",
        "is 0:",
        empty,
        ".\nDiscarding the sequence(s).\n"
      ))
      x <- x[-empty]
    }
    ctb <- .codonTable(x)
    ctb <- ctb[,order(colnames(ctb))] # sort codons alphabetically
    if (is.integer(ctb)) # in case there is only one sequence
      ctb <- rbind(ctb)
    namesx <- names(x)
    namesx[is.na(namesx)] <- "NA"
    new(
      "codonTable",
      ID = if (!is.null(namesx)) namesx else character(),
      counts = ctb,
      len = unname(rowSums(ctb)),
      KO = stringr::str_extract(namesx,"K\\d{5}"),
      COG = stringr::str_extract(namesx,"([KCN]|TW)OG\\d{5}")
    )
  }
)

#' @rdname codonTable-class
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
    
    rownamesx <- rownames(x)
    new(
      "codonTable",
      ID = if (!is.null(rownamesx)) rownamesx else character(),
      counts = rbind(x[,sort(colnames(x))]), # sort codons alphabet.
      len = rowSums(as.matrix(x), na.rm = TRUE),
      KO = stringr::str_extract(rownamesx,"K\\d{5}"),
      COG = stringr::str_extract(rownamesx,"([KCN]|TW)OG\\d{5}")
    )
  }
)

#' @rdname codonTable-class
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
    cat("codonTable instance with codon counts from", ns, "sequences.\n")
    df <- as.data.frame(object@counts)
    if (length(object@KO) != 0 & any(!is.na(object@KO)))
      df <- cbind(KO = object@KO, df, stringsAsFactors = FALSE)
    if (length(object@COG) != 0 & any(!is.na(object@COG)))
      df <- cbind(COG = object@COG, df, stringsAsFactors = FALSE)
    df <- cbind(length = object@len, df)
    if (length(object@ID) != 0 & any(!is.na(object@ID)))
      df <- cbind(ID = object@ID, df, stringsAsFactors = FALSE)
    if (ns < 15) {
      print(df)
    } else {
      out <- rbind(head(df), "..." = " ", tail(df))
      print(out)
    }
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
#' Get KO annotations of sequences
#' in \code{codonTable} object.
#'
#' @inheritParams codonCounts
#'
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
#' @describeIn codonTable Set KO annotations
#' for \code{codonTable} object.
#'
#' @inheritParams codonCounts
#' @param ann A character vector of sequence annotations,
#'    must be of length equal to \code{length(object)}.
#'
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
#' Get COG annotations of sequences
#' in \code{codonTable} object.
#'
#' @inheritParams codonCounts
#'
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

#' @describeIn codonTable Set COG annotations
#' for \code{codonTable} object.
#'
#' @inheritParams setKO
#'
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

#' @inheritParams base::Extract
#' @rdname subset
#' @export
setMethod(
  f="[",
  signature = "codonTable",
  definition = function(x,i){
    new("codonTable",
        ID = x@ID[i],
        counts = rbind(x@counts[i, ]),
        len = x@len[i],
        KO = x@KO[i],
        COG = x@COG[i])
  })

#' @rdname subset
#' @export
setMethod(
  f="[[",
  signature = "codonTable",
  definition = function(x,i){
    new("codonTable",
        ID = x@ID[i],
        counts = rbind(x@counts[i, ]),
        len = x@len[i],
        KO = x@KO[i],
        COG = x@COG[i])
  })

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
#' @return subsets of \code{codonTable} object, keeping in each slot
#' only those elements that meet the criteria in \code{subset}, if specified.
#'
#' @examples
#' # create codonTable
#' mat <- matrix(sample(1:10, 610, replace = TRUE), nrow = 10)
#' cT <- codonTable(mat) # produces informative warning
#' cT
#' cT[1]
#' cT[[1]]
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
