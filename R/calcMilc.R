# optimized MILC calculator
# the only required parameter is the codon usage table set,
# and MILC will be calculated against the average CU of that set.
# one optional parameter:
# subsets - (named) list of
#         - logical vectors of the same length as nrow(s), or
#         - additional CU table set(s) that will be processed
# the function returns a data frame of the same row length as nrow(s)
# and ncol equal to the number of subsets plus one (self)

calcMilc<-function(s, subsets=list()) {

    # test parameters
    if (length(subsets$ribosomal)>1) if (all(!subsets$ribosomal)) warning("Subset is empty! Please provide a valid subset!")
    if(! inherits(s, "codonTable")) stop("First argument must be a codon usage table!")
    nseq <- nrow(s)

    if(!is.list(subsets)) stop(paste("subsets must be a (named) list of logical vectors, each of length", nseq,
                                     "or codonTable objects (of any length)"))
    if(length(subsets) != 0) {
        ok <- sapply(subsets, function(x) {
            all(is.vector(x, mode = "logical"), length(x) == nseq) |
                all(inherits(x, "codonTable"), nrow(x) > 0)
        })
        stopifnot(ok)
        nam <- names(subsets)
        nsubs <- length(subsets)
        if(is.null(nam)) {
            nam <- paste("subset", 1:nsubs, sep = ".")
        } else {
            nam[nam == ""] <- paste("subset", (1:nsubs)[nam == ""], sep = ".")
        }
        names(subsets) <- make.names(nam, unique = TRUE)
    }

    # add a dummy self selection to subsets
    self_set <- rep(TRUE, nseq)
    subsets <- c(list(self = c(self_set)), subsets)

    # strip all unneeded info and convert to matrix
    o<-as.matrix(s[,nostops])

    # preprocess all subsets
    gc_list <- sapply(subsets, function(y) {
        if(is.vector(y, mode = "logical")) {
            sel <- o[y,] # this should never give error, because we tested for equal length
        } else {
            sel <- as.matrix(y[,nostops])
        }
        sel_sum <- colSums(sel) # add counts per codon
        gc <- byaa(sel_sum) # and normalize synonymous codons to sum = 1
        gc
    }, simplify = FALSE)

    # make fc values from original table
    fc<-byaas(o)

    # and correction factor
    corr1 <- t(t(byrcs(fc)) * as.vector(acnt-1))
    l <- s$len
    cf <- rowSums(corr1)/l - 0.5

    # now loop through all the gc's and calculate distance for each gene in the original set
    milcs <- sapply(gc_list, function(gc) {

        ma <- 2 * o[,order(ctab$aa)] * log(t(t(fc)/gc))
        MILC <- rowSums(ma, na.rm = TRUE)/l - cf
    })
    if(sum(s$len==0)) {s[s$len==0,"problem"]<-TRUE
    warning("Sequences of length 0 exist")
    }
    return(cbind(s, milcs))
}
