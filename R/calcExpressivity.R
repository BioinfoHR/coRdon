#' Calculation of CU-based expression level predictors.
#'
#' Calculates statistics that predict relative gene expression levels,
#' based on different synonymous codon usage (CU).
#'
#' @inheritParams calcCU
#' @param perc Top percent of sequences to be used as reference set in calculation
#'   of GCB statistic.
#'
#' @import data.table
#' @export
calcExpressivity <- function(cdt,
                             method,
                             subsets=list(),
                             seed = NULL,
                             ribosomal = FALSE,
                             gencode = "1",
                             alt.init = TRUE,
                             stop.rm = FALSE,
                             perc = 0.02) {

    if (!inherits(cdt, "codonTable"))
        stop("First argument must be a codon usage table!")

    if (!is.data.table(cdt))
        cdt <- as.data.table(cdt)

    codontab <- getGenCode(gencode, alt.init)
    ctab <- codontab$ctab
    stops <- codontab$stops
    nostops <- codontab$nostops
    len <- cdt$len.stop
    if (stop.rm |
        all(stops %in% colnames(cdt)) == FALSE |
        method %in% c("CAI")) {
        ctab <- ctab[ctab$codon %in% nostops, ]
        len <- cdt$len
    }
    cl <-
        lapply(levels(droplevels(ctab$AA)), function(x)
            which(ctab$AA == x))
    counts <- cdt[, ctab$codon, with = FALSE]
    csums <-
        sapply(cl, function(x)
            counts[, Reduce('+', .SD), .SDcols = x])
    freqs <- sapply(seq_along(cl), function(x)
        counts[, .SD / csums[, x], .SDcols = cl[[x]]])
    fc <-
        setDT(unlist(freqs, recursive = FALSE), check.names = TRUE)
    setcolorder(fc, ctab$codon)

    if (method %in% c("CAI","Fop")) {
        gc <- expectCU(cdt, ctab, subsets, self = FALSE, ribosomal)
        gcmax <-
            apply(gc, 2, function(x) # freq of the most frequent codons by aa
                by(x, droplevels(ctab$AA), max))
        x <- as.integer(droplevels(ctab$AA))
        gcm <- sapply(colnames(gc), function(j) # freq of the most frequent codons
            gc[, j] <- gcmax[x, j])

        if (method == "CAI"){
            cai <- lapply(colnames(gcm), function(y)
                exp(rowSums(counts * t(log(
                    t(fc) / gcm[, y] # relative adaptivenes of codons
                )), na.rm = TRUE) / len))
            cdt[, paste("CAI", colnames(gcm), sep = "_") := cai]

        } else if (method == "Fop"){
            fop <- lapply(colnames(gcm), function(y){
                ra <- t(fc) / gcm[, y]
                cnt <- as.matrix(counts)
                cnt[which(t(ra) < 0.9, arr.ind = T)] <- NA # as implemented in INCA
                rowSums(cnt, na.rm = T) / len
            })
            cdt[, paste("Fop", colnames(gcm), sep = "_") := fop]
        }

    } else if (method == "GCB"){

        if (length(seed) == 0 & ribosomal == FALSE)
            stop("Seed is not specified!")

        # seed <- sample(c(TRUE, FALSE), nrow(cdt), replace = TRUE)
        gcb_prev <- numeric(nrow(cdt))
        gc <- expectCU(cdt,
                       ctab,
                       list(seed),
                       self = FALSE,
                       ribosomal)
        iter <- 0
        repeat {
            cb <- log(gc / colMeans(fc, na.rm = T))
            cb[gc == 0] <- -5
            gcb <- rowSums(t(t(counts) * as.numeric(cb)), na.rm = T) / len
            diff <- all(gcb == gcb_prev)
            if (diff | iter > 6)
                break
            else {
                iter <- iter + 1
                gcb_prev <- gcb
                top <- order(gcb, decreasing = T)[1:(perc * length(gcb))]
                gc <- expectCU(cdt,
                               ctab,
                               list(seed = counts[top,]),
                               self = FALSE,
                               ribosomal)
            }
        }
        cdt[, GCB := gcb]

    } else if (method == "E"){
        calcCU(cdt, "B", subsets, gencode, stop.rm, alt.init, self = TRUE, ribosomal)
        es <- lapply(names(subsets), function(y) cdt$self / cdt[[y]])
        cdt[, paste("E", names(subsets), sep = "_") := es]

    } else if (method == "MELP"){
        calcCU(cdt, "MILC", subsets, gencode, stop.rm, alt.init, self = TRUE, ribosomal)
        melp <- lapply(names(subsets), function(y) cdt$self / cdt[[y]])
        cdt[, paste("MELP", names(subsets), sep = "_") := melp]
    }
}
