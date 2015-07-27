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
