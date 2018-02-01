#'
#' Calculate effective number of codons.
#'
#' Calculates effective number of codons (ENC) given set(s) of CU homozygosity
#' values.
#'
#' @param fa a vector containing homozygosity values for all 20 amino acids, or
#'   a matrix with homozygosity values for aminno acids in columns. If \code{fa} is
#'   a matrix, ENC is calculated for every row.
#' @param deg
#'
effNc <- function(fa, deg) {
    red <- unique(deg)[unique(deg) != 1]
    favg <- sapply(red, function(x){
        cols <- which(deg == x)
        avg <- rowSums(fa[,.SD, .SDcols=cols], na.rm = T) /
            rowSums(fa[, .SD, .SDcols=cols] != 0, na.rm = T)
        # if any redundancy class but F3 is absent
        if (x != 3)
            avg <- replace(avg, which(is.na(avg)), 1/x)
        length(cols) / avg
    })
    colnames(favg) <- red
    # if F3 is absent, average F2 and F4
    if (any(is.na(favg[,"3"]))) {
        n <- which(is.na(favg[,"3"]))
        if (all(favg[n, c("2","4")] > 0)) # ENC will be NA otherwise
        favg[n,"3"] <- sapply(n, function(i)
            (favg[i, "2"]/sum(deg == 2) + favg[i, "4"]/sum(deg == 4)) / 2)
    }
    enc <- sum(deg==1)+rowSums(favg)
    replace(enc, which(enc>61), 61)
}
