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
