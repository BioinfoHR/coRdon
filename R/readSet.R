#' Read Set
#'
#' Reads sets of fasta sequences stored in \code{folder}.
#'
#' @return Returns codon usage table for all files.
#'
#' @import data.table
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom purrr map
#' @importFrom dplyr progress_estimated
#' @export
#'
readSet <- function(folder = ".",
                    KOs = c(),
                    gencode = "1",
                    zipped = FALSE) {
  if (length(KOs) == 0)
    pattern <- "(*.fasta|*.FASTA)$"
  else
    pattern <- paste(KOs, sep = "|")

  if (zipped) {
    where <- tempdir()
    unzip(folder, exdir = where)
    folder <- where
  }
  files <- dir(folder, pattern = pattern)

  combined <- map(files, function(x) {
    aset <- readDNAStringSet(paste(folder, x, sep = "/"))
    as.character(aset)
  })
  names(combined) <- files

  if (zipped)
    unlink(where, recursive = TRUE)

  DNAStringSet(unlist(combined))
}
