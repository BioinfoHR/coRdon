readSet <- function(folder = ".", KOs = c(), zipped = FALSE) {

  if(length(KOs) == 0)
    pattern <- "*.fasta"
  else
    pattern <- paste(KOs, sep="|")

  if(zipped) {
    where <- tempdir()
    unzip(folder, exdir = where)
    folder <- where
  }
  files <- dir(folder, pattern=pattern)

  combined <- llply(files, function(x) {
      aset <- readDNAStringSet(paste(folder, x, sep="/"))
      as.character(aset)
    }, .progress = "text")
  names(combined) <- files

  if(zipped) unlink(where, recursive = TRUE)

  aset <- DNAStringSet(unlist(combined))

  ctable <- oligonucleotideFrequency(aset, width = 3, step = 3)

#  KO <- str_replace_all(names(aset), ".*(K\\d{5}).*", "\\1")
#  COG <- str_replace_all(names(aset), ".*(([KCN]|TW)OG\\d{5}).*", "\\1")

  ID <- names(aset)
  KO <- str_extract(ID, "K\\d{5}")
  COG <- str_extract(ID, "([KCN]|TW)OG\\d{5}")
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
