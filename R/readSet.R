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
#' @param prepend.filenames Logical, whether to prepend filename(s) 
#'    to names in \code{DNAStringSet} object. Default is \code{FALSE}.
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
#' @importFrom dplyr progress_estimated
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom purrr map
#'
#' @export
readSet <- function(
    folder = character(), file = character(),
    KOs = c(), zipped = FALSE, prepend.filenames = FALSE) {
    
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
        names <- dir(folder, pattern = pattern)
        files <- dir(folder, pattern = pattern, full.names = TRUE)
    } else if (length(file)!=0) {
        if (zipped) {
            where <- tempdir()
            unzip(file, exdir = where)
            folder <- where
            names <- dir(folder, pattern = pattern)
            files <- dir(folder, pattern = pattern, full.names = TRUE)
        } else {
            names <- file
            files <- file
        }
    }
    
    pb <- progress_estimated(length(files), min_time = 1)
    
    combined <- map(files, function(x) {
        pb$tick()$print()
        aset <- readDNAStringSet(x)
        as.character(aset)
    })
    
    if (prepend.filenames) {
        names(combined) <- names
    }
    
    if (zipped)
        unlink(where, recursive = TRUE)
    
    DNASS <- DNAStringSet(unlist(combined))
    return(DNASS)
}
