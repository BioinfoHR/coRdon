#' Calculation of CU statistics.
#'
#' Calculates codon usage (CU) statistics for given codonTable object.
#'
#' @param cdt A codonTable object. Rows represent indivisual sequences. Must
#'   contain the following columns:
#'   \itemize{
#'     \item len, sequence length in codons, stop codons excluded
#'     \item len.stop, sequence length in codons, including stop codons
#'   }
#'   and a single column with counts of every codon.
#' @param method A character string indicating which CU statistic will be
#'   calculated, one of the following: "MILC", "B", "MCB", "ENC", "ENCp", "SCUO".
#' @param subsets A (named) list of logical vectors, the length of each equal
#'   to the number of sequences in the set, or character vectors (of any length)
#'   containing KEGG/eggNOG annotations, or codonTable objects (of any length).
#' @param gencode A character that uniquely identifies the genetic code to be used
#'   for calculation of CU statistic. Should be one of the values in the \code{id}
#'   or \code{name2} columns of \code{GENETIC_CODE_TABLE}.
#' @param alt.init Logical, whether to include alternative initiation codons in
#'   the genetic code.
#' @param stop.rm Logical, whether to include stop codons (if they are present)
#'   in CU calculation. Default is \code{FALSE}, note however that if method is
#'   ENC or ENCp, stop codons are not used by definition.
#' @param self Logical, if \code{TRUE} (default), CU statistic is also calculated
#'   against the average CU of the entire sequence set. Just as \code{subsets},
#'   this argument is also intended for CU statistics that account for background
#'   nucleotide composition.
#' @param ribosomal Logical, if \code{TRUE}, CU statistic is also calculated
#'   against the average CU of the ribosomal genes in the sequence set. This
#'   argument is also intended for CU statistics that account for background
#'   nucleotide composition.
#'
#' @details \code{subsets} should only be supplied if CU statistic that accounts
#'   for background nucleotide composition is calculated, i.e if \code{method} is
#'   one of the following: "MILC", "B", "MCB", "ENCp". Otherwise \code{subsets}
#'   parameter is ignored.
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
                   gencode = "1",
                   alt.init = TRUE,
                   stop.rm = FALSE,
                   self = TRUE,
                   ribosomal = FALSE) {

    if (!inherits(cdt, "codonTable"))
        stop("First argument must be a codon usage table!")

    if (!is.data.table(cdt))
        cdt <- as.data.table(cdt) # setDT makes a shallow copy within fun env

    codontab <- getGenCode(gencode, alt.init)
    ctab <- codontab$ctab
    stops <- codontab$stops
    nostops <- codontab$nostops
    len <- cdt$len.stop

    if (stop.rm |
        all(stops %in% colnames(cdt)) == FALSE |
        method %in% c("ENC", "ENCp")) {
        ctab <- ctab[ctab$codon %in% nostops, ]
        len <- cdt$len
    }

    cl <- lapply(levels(droplevels(ctab$AA)), function(x) which(ctab$AA==x))
    deg <- sapply(cl, length)
    counts <- cdt[, ctab$codon, with=FALSE]
    # codon counts summed by AA
    csums <- sapply(cl, function(x) counts[, Reduce('+',.SD), .SDcols = x])
    # normalized codon frequencies
    freqs <- sapply(seq_along(cl), function(x)
        counts[,.SD/csums[,x], .SDcols = cl[[x]]])
    fc <- setDT(unlist(freqs, recursive = FALSE), check.names = TRUE)
    setcolorder(fc, ctab$codon)

    if (method %in% c("MILC","B","MCB","ENCp")) {
        # expected codon frequencies
        gc <- expectCU(cdt, ctab, subsets, self, ribosomal)

        if (method == "MILC") {
            cor <- as.numeric(((csums>0) %*% (deg-1))/len - 0.5)
            milc <- lapply(colnames(gc), function(x) {
                rowSums(2 * counts * log(t(t(fc) / gc[, x])), na.rm = TRUE) / len - cor
            })
            cdt[, colnames(gc) := milc]

        } else if (method == "B") {
            pa <- csums / len
            b <- lapply(colnames(gc), function(x) {
                dt <- as.data.table(t(abs(t(fc) - gc[, x])))
                ba <- sapply(cl, function(x)
                    dt[, Reduce('+', .SD), .SDcols = x])
                rowSums(pa * ba, na.rm = TRUE)
            })
            cdt[, colnames(gc) := b]

        } else if (method == "MCB") {
            A <- csums[, deg>1] > 0
            mcb <- lapply(colnames(gc), function(x) {
                dt <- as.data.table(t((t(fc) - gc[, x])^2 / gc[, x]))
                dt[which(is.na(dt), arr.ind = T)] <- 0 # gc > 0
                dt[counts<=0] <- 0 # counts > 0
                ba <- sapply(cl, function(y)
                    dt[, Reduce('+', .SD), .SDcols = y])
                rowSums(ba[, deg>1] * log10(csums[, deg>1]), na.rm = T) / rowSums(A)
            })
            cdt[, colnames(gc) := mcb]


        } else if (method == "ENCp"){
            encp <- lapply(colnames(gc), function(x) {
                # chi squared
                dt <- as.data.table(t((t(fc) - gc[, x])^2 / gc[, x]))
                if (any(gc==0, na.rm = T)) {
                    infs <- unique(which(dt == Inf, arr.ind = T)[,2])
                    dt[, (infs) := 0]
                }
                dt[which(is.na(dt), arr.ind = T)] <- 0
                chisum <- sapply(cl, function(x)
                    dt[, Reduce('+', .SD), .SDcols = x])
                chisq <- csums * chisum
                # homozygosity
                fa <- as.data.table(t((t(chisq + csums) - deg) /
                                          (deg * t(csums - 1))))
                fa[csums<5] <- NA
                effNc(fa, deg)
            })
            cdt[, colnames(gc) := encp]
        }

    } else if (method == "ENC") {
        pi <- sapply(freqs, function(x) rowSums(x^2, na.rm = TRUE))
        # remove AA with count <= 1
        csums[which(csums <= 1, arr.ind = T)] <- NA
        fa <- as.data.table((csums * pi - 1) / (csums - 1))
        enc <- effNc(fa, deg)
        cdt[, ENC := enc]

    } else if (method == "SCUO") {
        Ha <- sapply(freqs, function(x) rowSums(-x*log10(x), na.rm = TRUE))
        Hmax <- log10(sapply(cl, length))
        Oa <- t((Hmax - t(Ha)) / Hmax)
        Fa <- csums / rowSums(csums[,deg>1])
        scuo <- rowSums(Oa * Fa, na.rm = TRUE)
        cdt[, SCUO := scuo]
    }

    return(cdt)
}
