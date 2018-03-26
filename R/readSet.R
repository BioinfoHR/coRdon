#' Read set of sequences.
#'
#' Reads set of fasta files stored in \code{folder}.
#'
#' @param folder Path to directory containing .fasta files.
#' @param KOs An optional character vector of sequence annotations (e.g. KO)
#'    contained in the names of fasta files to be selectively read.
#' @param zipped Logical, whether \code{folder} is zipped. Default is
#'    \code{FALSE}.
#'
#' @return Returns a \code{DNAStringSet} object.
#'
#' @examples
#' exampledir <- system.file("extdata", package = "coRdon")
#' list.files(exampledir)
#' readSet(exampledir)
#' readSet(exampledir, KOs = "K02438")
#'
#' @import data.table
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom purrr map
#'
#' @export
readSet <- function(folder = ".",
                    KOs = c(),
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
