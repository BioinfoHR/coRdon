library(Biostrings)
library(Biobase)
library(plyr)
library(stringr)
#library(Rcpp)
pldensity<-function(a,groupnames){
    #a is aggregated data.frame BY NAMES (H1,H2,...,L1,...)
    a<-a[!is.na(a$x),]
    l <- list()
    y <- list()
    for (i in 1:length(groupnames)){
        l[[i]] <- grep(groupnames[i],a$Group.1)
        if (length(l[[i]])) y[[i]] <- density(a$x[l[[i]]])
        #lines(y,col=i)
    }
    x_lim <- range(sapply(y, function(x){ range(x$x)}))
    y_lim <- range(sapply(y, function(x){ range(x$y)}))

    plot(y[[1]],ylim=y_lim,xlim=x_lim,main = "Density by groups")
    for (i in 2:length(y)){
        lines(y[[i]],col=i)
    }
}


#################################################################




min_length <- 30
percentile_top <- 0.90
perc <- c(0.95, 0.90, 0.85, 0.70, 0.50)
#KEGG_maps_folder <- "KEGG_maps"
#location_of_ko <- "C:/Users/mfabijanic/Dropbox/Dusko/maja/ko.Robj"
#location_of_brite <- "C:/Users/mfabijanic/Dropbox/Dusko/maja/brite.Robj"
# load KO and BRITE ontologies (use separate script to prepare)
#load(location_of_ko)
#load(location_of_brite)

reduce.contable <- function(ctb, column = "B") {
    values <- unique(ko[,column])
    tt <- lapply(values, function(x){
        KOs <- ko[ko[,column] == x, "KO", drop=TRUE]
        KOs <- unique(KOs)
        counts <- subset(ctb, KO %in% KOs , select = -KO)
        colSums(counts, na.rm = TRUE)
    })
    contable <- do.call(rbind, tt)
    rownames(contable) <- values
    as.data.frame(contable)
}

make.stats <- function(contable) {

    rows <- names(contable)
    top_rows <- rows[grep("top", rows)]

    all <- contable$all
    all.sum <- sum(all)

    by.top <- sapply(top_rows, function(row) {

        top <- contable[,row]
        top.sum <- sum(top)

        ct <- data.frame(all = all, top)
        names(ct) <- c("all", "cnt")
        rownames(ct) <- rownames(contable)
        sc <- top.sum / all.sum
        scaled_top <- top + 1
        scaled_all <- all * sc + 1
        ct$enrich <- (scaled_top - scaled_all) / scaled_all * 100
        ct$M <- log2(scaled_top) - log2(scaled_all)
        ct$A <- (log2(scaled_all) + log2(scaled_top)) / 2
        pvals =
            apply(ct[,c("all", "cnt")], 1, function(x) {
                b = binom.test(x[2], top.sum, x[1]/all.sum)
                b$p.value
            })
        ct$pvals = pvals
        ct$padj = p.adjust(pvals, method = "BH")
        ct$all <- NULL
        ct

    }, simplify = FALSE)

    bound <- do.call(cbind, by.top)
    data.frame(all = all, bound)

}



make.contable <- function(csm, KEGG) {
    if (!KEGG) csm$KO <- as.factor(as.character(csm$ID))
    csm$KO <- as.factor(ifelse(!is.na(as.character(csm$COG)),as.character(csm$COG),as.character(csm$KO)))
    all <- as.vector(table(csm$KO))
    top <- sapply(perc, function(x) {
        as.vector(table(csm[csm$melp >= quantile(csm$melp, x), "KO"]))
    })
    top <- as.data.frame(top)
    names(top) <- make.names(paste("top", perc, sep="_"))

    data.frame(all = all, top, row.names = levels(csm$KO), KO = levels(csm$KO))
}
doall_csm <- function(csm,levels_KEGG=c()){
    csm <- csm[csm$len > min_length,]
    enrichment <-list()
    KEGG <- T
    if (all(is.na(csm$KO))) KEGG <- F
    if(KEGG)
        if(!(exists("ko")&&exists("brite"))) {
            KEGG <- F
            print("ko.Robj and brite.Robj not found")
        }
    if(!KEGG) csm$KO <- csm$ID
    # three levels of ontology exist: B (very broad), C (pathways), and KO (orthology)
    # we'll summarize by each level separately

    # make a contingency table for KO level, this is the base for other two
    contable_ko <- make.contable(csm, KEGG)
    name <- unique(csm$name)

    if (KEGG){

        csm$KO <- as.factor(csm$KO)

        if ("B" %in% levels_KEGG){
            contable_b <- reduce.contable(contable_ko, "B")
            contable_b <- make.stats(contable_b)
            contable_b$B <- rownames(contable_b)
            #    write.table(contable_b, file=paste(st,"contable_b.txt", sep="_"), row.names = FALSE, sep = "\t")
            enrichment[[name]][["B"]] <- contable_b
            print(paste("Enrichment calculated for level B for name ",name,sep=""))
        }
        if ("C" %in% levels_KEGG){

            contable_c <- reduce.contable(contable_ko, "C")
            contable_c <- make.stats(contable_c)
            contable_c$C <- rownames(contable_c)
            contable_cc <- merge(contable_c, brite, by.x="C", by.y="C", all.x=TRUE)
            contable_cc$name <- name
            #    contable_cc <- contable_cc[order(contable_cc$padj),]
            #    write.table(contable_cc, file=paste(st,"contable_c.txt", sep="_"), row.names = FALSE, sep = "\t")
            enrichment[[name]][["C"]] <- contable_cc
            print(paste("Enrichment calculated for level C for name ",name,sep=""))
        }
    }
    contable_koc <- make.stats(contable_ko)
    contable_koc$KO <- rownames(contable_koc)
    if (KEGG){
        contable_kocc <- merge(contable_koc, unique(ko[,c("KO", "desc")]), by.x="KO", by.y="KO", all.x=TRUE)
        kos <- unique(csm$KO)
    }else{
        contable_kocc <- contable_koc
        kos <- unique(csm$ID)
    }
    contable_kocc$median <- 0
    if (!KEGG)
    {
        contable_kocc$desc <- name
        contable_kocc$name <- name
        for (i in 1:length(kos))
            contable_kocc[contable_kocc$KO==kos[i],"median"] <- median(csm[csm$ID==kos[i],"melp"])
    }
    if (KEGG)
    {
        contable_kocc$name <- name
        medianis2 <- aggregate.data.frame(csm$melp,by=list(csm$KO),FUN=median)
        contable_kocc$median <- 0
        for(i in medianis2$Group.1)
        {
            contable_kocc[contable_kocc$KO==i,]$median <- medianis2[medianis2$Group.1==i,]$x
        }
    }
    enrichment[[name]][["KO"]] <- contable_kocc
    #    contable_kocc <- contable_kocc[order(contable_kocc$padj),]
    #    write.table(contable_kocc, file=paste(st,"contable_ko.txt", sep="_"), row.names = FALSE, sep = "\t")
    if(KEGG) enrichment[[name]][["KO"]]<-enrichment[[name]][["KO"]][,c(2:32,1,33:35)]
    enrichment[[name]][["KO"]] <- enrichment[[name]][["KO"]][!is.na(enrichment[[name]][["KO"]]$median),]
    enrichment[[name]][["KO"]] <- enrichment[[name]][["KO"]][!duplicated(enrichment[[name]][["KO"]]$KO),]
    if ("C"%in% levels_KEGG)
    {
        #enrichment[[name]][["C"]] <- enrichment[[name]][["C"]][!is.na(enrichment[[name]][["C"]]$median),]
        enrichment[[name]][["C"]] <- enrichment[[name]][["C"]][!duplicated(enrichment[[name]][["C"]]$C),]
        enrichment[[name]][["C"]]$name <- name
    }
    if ("B"%in% levels_KEGG) {
        #enrichment[[name]][["B"]] <- enrichment[[name]][["B"]][!is.na(enrichment[[name]][["B"]]$median),]
        enrichment[[name]][["B"]] <- enrichment[[name]][["B"]][!duplicated(enrichment[[name]][["B"]]$B),]
        enrichment[[name]][["B"]]$name <- name

    }
    return(enrichment[[name]])
}
