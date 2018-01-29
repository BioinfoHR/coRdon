#' Calculation of expected codon frequencies.
#'
#' Calculates expected codon frequencies based on CU distribution for each
#' \code{subsets} of sequences in \code{cdt}.
#'
#' @inheritParams calcCU
#'
#' @return A 61 x \code{length(subsets)} matrix with alphabetically ordered codons'
#'   frequencies (rows) for every subset (column).
#'
#' @import data.table
#'
expectCU <- function(cdt, ctab, subsets, self, ribosomal){

    ns <- nrow(cdt)

    if (!is.list(subsets))
        stop(
            paste(
                "Subsets must be either a (named) list of logical vectors, each
                of length", ns,
                "or character vectors containing KEGG/eggNOG annotations,
                or codonTable objects (of any length)."
            )
        )
    if (length(subsets) != 0) {
        ok <- sapply(subsets, function(x) {
            all(is.vector(x, mode = "logical"), length(x) == ns) |
                all(grep("K\\d{5}", x)) |
                all(grep("([KCN]|TW)OG\\d{5}", x)) |
                all(inherits(x, "codonTable"), nrow(x) > 0)
        })
        stopifnot(ok)
        nam <- names(subsets)
        nsubs <- length(subsets)
        if (is.null(nam)) {
            nam <- paste("subset", 1:nsubs, sep = ".")
        } else {
            nam[nam == ""] <- paste("subset", (1:nsubs)[nam == ""], sep = ".")
        }
        names(subsets) <- make.names(nam, unique = TRUE)
    }

    if (ribosomal == TRUE) {
        subsets <- c(list(ribosomal = cdt$KO %in% RPKOs), subsets)
    }
    if (self == TRUE) {
        self_set <- rep.int(TRUE, ns)
        subsets <- c(list(self = c(self_set)), subsets)
    }

    sapply(subsets, function(x) {
        if (is.vector(x, mode = "logical")) {
            x <- cdt[x, ctab$codon, with=FALSE]
        } else if ( all(grep("K\\d{5}", x))) {
            x <- cdt[KO %in% x, ctab$codons, with=FALSE]
        } else if (all(grep("([KCN]|TW)OG\\d{5}", x))) {
            x <- cdt[COG %in% x, ctab$codon, with=FALSE]
        } else {
            x <- as.data.table(x)[, ctab$codon, with=FALSE]
        }
        x <- x[, lapply(.SD, sum), .SDcols = ctab$codon]
        melt(x,
             measure.vars = ctab$codon,
             variable.name = "codon",
             value.name = "c"
        )[, aa := ctab$aa][, g := c/sum(c), by = aa][, g]
    })
}
