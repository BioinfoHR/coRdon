#' @import KEGGREST
#' @import data.table

paths <- names(keggList("pathway"))
paths <- regmatches(paths, regexpr("[[:alpha:]]{2,4}\\d{5}", paths)) # 523 paths
pnames <- unname(keggList("pathway"))
kop <- lapply(1:length(paths), function(x){
    KO <- unname(keggLink("ko", paths[x]))
    KO <- regmatches(KO, regexpr("K\\d{5}", KO))
    if (length(KO) == 0) KO <- NA
    data.table(CATEGORY = paths[x], ANN = KO, description = pnames[x])
})
KO_PATHWAYS <- rbindlist(kop)

mods <- names(keggList("module"))
mods <- regmatches(mods, regexpr("M\\d{5}", mods)) # 788 modules
mnames <- unname(keggList("module"))
kom <- lapply(1:length(mods), function(x){
    KO <- unname(keggLink("ko", mods[x]))
    KO <- regmatches(KO, regexpr("K\\d{5}", KO))
    if (length(KO) == 0) KO <- NA
    data.table(CATEGORY = mods[x], ANN = KO, description = mnames[x])
})
KO_MODULES <- rbindlist(kom)

COGs <- read.table("ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab",
                   sep = "\t",
                   header = FALSE,
                   comment.char = "#",
                   quote="",
                   fill=FALSE,
                   stringsAsFactors = FALSE)
COGs <- as.data.table(COGs)
setnames(COGs, c("V1","V2","V3"), c("ANN","CAT","description"))
mnchar <- max(sapply(unique(COGs$CAT), function(x) nchar(x), USE.NAMES = FALSE))
COGs[, paste0("V", 1:mnchar) := tstrsplit(CAT, "")]
COGs <- melt.data.table(COGs, measure.vars = paste0("V", 1:mnchar),
                        value.name = "CATEGORY", na.rm = TRUE)
COGs[, `:=`(CAT=NULL, variable=NULL)]
setcolorder(COGs, c("CATEGORY","ANN","description"))
setkey(COGs, ANN)

# devtools::use_data(KO_PATHWAYS, KO_MODULES, COGs, internal = TRUE, overwrite = TRUE)

dnaLD94 <- readSet(file = "http://hex.bioinfo.hr/~mfabijanic/LD94.fasta")
LD94 <- codonTable(dnaLD94[1001:2000])
dnaHD59 <- readSet(file = "http://hex.bioinfo.hr/~mfabijanic/HD59.fasta")
HD59 <- codonTable(dnaHD59[1001:2000])

# devtools::use_data(HD59, LD94, overwrite = TRUE)

s <- getKO(HD59)
v <- as.numeric(MELP(HD59, ribosomal = TRUE))
ct <- crossTab(s, v)
HD59_KO <- enrichment(ct)
ctp <- reduceCrossTab(ct, "pathway")
HD59_PATHWAYS <- enrichment(ctp)

s <- getKO(LD94)
v <- as.numeric(MELP(LD94, ribosomal = TRUE))
ct <- crossTab(s, v)
LD94_KO <- enrichment(ct)
ctp <- reduceCrossTab(ct, "pathway")
LD94_PATHWAYS <- enrichment(ctp)

#devtools::use_data(HD59_KO, HD59_PATHWAYS, LD94_KO, LD94_PATHWAYS, overwrite = TRUE)
