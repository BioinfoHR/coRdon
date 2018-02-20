#' @include codonTable-class.R
#' @include genCode-class.R
#' @import data.table

.countsbyaa <- function(cTobject, gCobject)
{
    counts <- as.data.table(cTobject@counts)
    cl <- gCobject@cl
    rbind(sapply(cl, function(x) counts[, Reduce('+',.SD), .SDcols = x]))
}

.normFrequencies <- function(cTobject, gCobject)
{
    counts <- as.data.table(cTobject@counts)
    csums <- .countsbyaa(cTobject, gCobject)
    cl <- gCobject@cl
    freqs <- sapply(seq_along(cl), function(x)
        counts[,.SD/csums[,x], .SDcols = cl[[x]]])
    fc <- setDT(unlist(freqs, recursive = FALSE))
    setcolorder(fc, sort(gCobject@codons)) # sort codons alphabetically
    return(fc)
}

.normSetFrequencies <- function(cTobject, gCobject)
{
    csum <- as.data.table(rbind(colSums(cTobject@counts)))
    csumaa <- colSums(.countsbyaa(cTobject, gCobject))
    cl <- gCobject@cl
    out <- sapply(seq_along(cl), function(x)
        csum[, .SD/csumaa[x], .SDcols = gCobject@cl[[x]]])
    unlist(out)[sort(gCobject@codons)] # sort codons alphabetically
}

.effNc <- function(fa, deg) {
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
        if (all(favg[n, c("2","4")] > 0))
            favg[n,"3"] <- sapply(n, function(i)
                (favg[i, "2"]/sum(deg == 2) + favg[i, "4"]/sum(deg == 4)) / 2)
        else favg[,"3"] <- 1/sum(deg == 3)
    }
    enc <- sum(deg==1)+rowSums(favg)
    replace(enc, which(enc>61), 61)
}
