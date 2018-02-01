#' Get genetic code table
#'
#' @importFrom Biostrings getGeneticCode
#' @import data.table
#'
getGenCode <- function(gencode, alt.init=TRUE) {
    ctab <- as.data.table(getGeneticCode(gencode, as.data.frame = T),
                          keep.rownames = "codon")
    ctab[, AA := as.factor(ctab$AA)]
    if (alt.init) ctab[Start == "M", AA:="M"]
    ctab[, Start := NULL]
    list(
        ctab = ctab,
        stops = ctab[AA == "*", codon],
        nostops = ctab[AA != "*", codon]
    )
}
