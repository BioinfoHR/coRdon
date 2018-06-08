#' Read set of sequences.
#'
#' Reads a set of fasta files stored in \code{folder},
#' or a single fasta \code{file}.
#'
#' @param folder Path to directory containing .fasta files.
#' @param file Path to a single .fasta file, or zipped file
#'     (if latter, specify \code{ZIPPED = TRUE}).
#' @param KOs An optional character vector of sequence annotations (e.g. KO)
#'    contained in the names of fasta files to be selectively read.
#' @param zipped Logical, whether \code{folder} or \code{file} is zipped.
#'    Default is \code{FALSE}.
#'
#' @return Returns a \code{DNAStringSet} object.
#'
#' @examples
#' exampledir <- system.file("extdata", package = "coRdon")
#' files <- list.files(exampledir)
#' readSet(folder = exampledir)
#' readSet(folder = exampledir, KOs = "K02931")
#' pathtofile <- paste(exampledir, files[1], sep = "/")
#' readSet(file = pathtofile)
#'
#' @import data.table
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom purrr map
#'
#' @export
readSet <- function(folder = character(),
                    file = character(),
                    KOs = c(),
                    zipped = FALSE) {
  if (length(KOs) == 0)
    pattern <- "(*.fasta|*.FASTA)$"
  else
    pattern <- paste(KOs, sep = "|")

  if (length(folder)!=0) {

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

  } else if (length(file)!=0) {

      if (zipped) {
          where <- tempdir()
          unzip(file, exdir = where)
          folder <- where

          files <- dir(folder, pattern = pattern)

          combined <- map(files, function(x) {
              aset <- readDNAStringSet(paste(folder, x, sep = "/"))
              as.character(aset)
          })
          names(combined) <- files

          unlink(where, recursive = TRUE)

          DNAStringSet(unlist(combined))

      } else {

          aset <- readDNAStringSet(file)
          combined <- list(as.character(aset))
          names(combined) <- regmatches(
              file, gregexpr("([^///]*.fasta|[^///]*.FASTA)$", file))[[1]]
          DNAStringSet(unlist(combined))

      }
  }
}
