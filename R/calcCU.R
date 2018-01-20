#' Calculation of CU statistics.
#'
#' Calculates codon usage (CU) statistics for given codonTable object.
#'
#' @param cdt A codonTable object. Rows represent indivisual sequences. Must
#'   contain the following columns:
#'   \itemize{
#'     \item ID, a unique sequence identifier
#'     \item len, sequence length in codons
#'   }
#'   and a single column with counts of every no-stop codon.
#' @param method A character string indicating which CU statistic will be
#'   calculated, one of the following: "MILC", "B", "MCB", "ENC", "ENCp", "SCUO".
#' @param subsets A (named) list of logical vectors, the length of each equal
#'   to the number of sequences in the set, or codonTable objects of any length.
#'   \code{subsets} should only be supplied if CU statistic that accounts for
#'   background nucleotide composition is calculated, i.e if \code{method} is
#'   one of the following: "MILC", "B", "MCB", "ENCp".
#' @param self Logical, if \code{TRUE} (default), CU statistic is also calculated
#'   against the average CU of the entire sequence set. Just as \code{subsets},
#'   this argument is also intended for CU statistics that account for background
#'   nucleotide composition.
#' @param ribosomal Logical, if \code{TRUE}, CU statistic is also calculated
#'   against the average CU of the ribosomal genes in the sequence set. This
#'   argument is also intended for CU statistics that account for background
#'   nucleotide composition.
#'
#' @return Returns an input codonTable object augmented with additional columns
#'   containing CU statistic values for each subset.
#'
#' @import data.table
#' @export
#'
calcCU <- function(cdt,
                   method,
                   subsets = list(),
                   self = TRUE,
                   ribosomal = FALSE) {

  len <- cdt$len
  counts <- cdt[, nostops, with=FALSE]
  # codon counts summed by aa
  csums <- sapply(cl, function(x) counts[, Reduce('+',.SD), .SDcols = x])
  # normalized codon frequencies
  freqs <- sapply(1:20, function(x) counts[,.SD/csums[,x], .SDcols = cl[[x]]])
  fc <- setDT(unlist(freqs, recursive = FALSE), check.names = TRUE)
  setcolorder(fc, ctab$codonstr)

  # expected codon frequencies
  gc <- expectCU(cdt, subsets, self)

  if (method == "MILC") {
      cor <- as.numeric(((csums>0) %*% (acnt-1))/len - 0.5)
      milc <- lapply(colnames(gc), function(x) {
          rowSums(2 * counts * log(t(t(fc) / gc[, x])), na.rm = TRUE) / len - cor
      })
      cdt[, colnames(gc) := milc]
  }

  return(cdt)
}
