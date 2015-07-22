library(Biostrings)
library(Biobase)
library(plyr)
library(stringr)

# helper constants

RPKOs = c("K02945", "K02967", "K02982", "K02986", "K02988", "K02990", "K02992", "K02994", "K02996",
          "K02946", "K02948" ,"K02950" ,"K02952" ,"K02954" ,"K02956" ,"K02959" ,"K02961" ,"K02963" ,
          "K02965" , "K02968" , "K02970" , "K02981" , "K02985" , "K02984" , "K02987" , "K02989" , "K02991",
          "K02993" , "K02995" , "K02997" , "K02947" , "K02949" , "K02951" , "K02953" , "K02955" , "K02958" ,
          "K02957" , "K02960" , "K02962" , "K02964" , "K02966" , "K02969" , "K02971" , "K02973" , "K02974" , "K02975" ,
          "K02976" , "K02978" , "K02977" , "K02979" , "K02980" , "K02983" , "K02998" , "K02863" , "K02886" ,
          "K02906" , "K02926" , "K02931" , "K02933" , "K02935" , "K02939" , "K02864" , "K02867" , "K02869" ,
          "K02871" , "K02874" , "K02876" , "K02878" , "K02879" , "K02881" , "K02884" , "K02887" , "K02888" ,
          "K02890" , "K02892" , "K02895" , "K02897" , "K02899" , "K02902" , "K02904" , "K02907" , "K02909" ,
          "K02911" , "K02913" , "K02914" , "K02916" , "K02919" , "K07590" , "K02925" , "K02930" , "K02932" ,
          "K02934" , "K02937" , "K02936" , "K02938" , "K02940" , "K02866" , "K02865" , "K02868" , "K02870" ,
          "K02873" , "K02872" , "K02875" , "K02877" , "K02880" , "K02883" , "K02882" , "K02885" , "K02889" ,
          "K02891" , "K02894" , "K02893" , "K02896" , "K02898" , "K02901" , "K02900" , "K02903" , "K02905" ,
          "K02908" , "K02910" , "K02912" , "K02915" , "K02918" , "K02917" , "K02920" , "K02922" , "K02921" ,
          "K02923" , "K02924" , "K02927" , "K02928" , "K02929" , "K02941" , "K02942" , "K02943" , "K02944" ,
          "K01977" , "K01980" , "K01985" , "K01979" , "K01982" , "K01981"
)

# constants

dummy<-"AAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTAATACTAGTATTCATCCTCGTCTTGATGCTGGTGTTTATTCTTGTTT"

codons<-strsplit(strbreak(dummy, width=3, exdent=0, collapse="|"), split="\\|")[[1]]
stops<-c("TAA", "TAG", "TGA")
nostops<-codons[!(codons %in% stops)]



ctab<-data.frame(
    codon=nostops,
    aa=as.factor(as.character(translate(DNAStringSet(nostops)))),
    codonstr=as.character(nostops), stringsAsFactors=F
)

acnt <- tapply(ctab$codon, ctab$aa, length)
ordaa <- ctab[order(ctab$aa),"aa"]

# functions

cnorm<-function(x) x/sum(x)
cnormm <- function(x) if(is.matrix(x)) x/rowSums(x) else x/x

rcm <- function(x) if(is.matrix(x)) as.integer(!is.nan(rowSums(x))) else as.integer(!is.nan(x))
byaa<-function(x) unlist(tapply(x,ctab$aa,cnorm))

byaas <- function(x) do.call(cbind, sapply(names(acnt), function(y) cnormm(x[,ctab$aa %in% y])))
byrcs <- function(x) sapply(names(acnt), function(y) rcm(x[,ordaa %in% y]))


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
    if (length(subsets)>1) if (all(!subsets)) warning("Subset is empty!")
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


readSet <- function(folder=., KOs=c(), zipped = FALSE){
    aset <- readDNAStringSet(folder)
    ctable <- oligonucleotideFrequency(aset, width = 3, step = 3)
    ID <- names(aset)
    KO <- str_extract(ID, "K\\d{5}")
    COG <- str_extract(ID, "([KCN]|TW)OG\\d{4}")
    len.stop <- rowSums(ctable[,codons])
    len <- rowSums(ctable[,nostops])
    problem <- rowSums(ctable) != len.stop

    ccc = data.frame(
        ID = ID,
        ctable,
        KO = KO,
        COG = COG,
        len.stop = len.stop,
        len = len,
        problem = problem
    )

    class(ccc) <- c(class(ccc), "codonTable")
    ccc
}
calculateMelp <- function(file,RPKOs){
    myset <- readSet(file)
    ribosomals <- myset$KO %in% RPKOs
    milc <- calcMilc(myset, list(ribosomal = ribosomals))
    milc$melp <- milc$self / milc$ribosomal
    milc$name <- file
    return(milc)
}
melp_all <- function(path){
    currwd <- getwd()
    setwd(path)
    sample <- list()
    sample <- lapply(list.files(path),function(x){calculateMelp(x,RPKOs)})
    all_files <- do.call("rbind",sample)
    setwd(currwd)
    rm(sample)
    all_files[!is.na(all_files$melp),]
    return(all_files)
}
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
    csm$KO <- csm$ID
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
