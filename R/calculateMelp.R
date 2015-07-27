calculateMelp <- function(file,RPKOs){
    myset <- readSet(file)
    ribosomals <- myset$KO %in% RPKOs
    milc <- calcMilc(myset, list(ribosomal = ribosomals))
    milc$melp <- milc$self / milc$ribosomal
    milc$name <- file
    return(milc)
}
