# functions

cnorm<-function(x) x/sum(x)
cnormm <- function(x) if(is.matrix(x)) x/rowSums(x) else x/x

rcm <- function(x) if(is.matrix(x)) as.integer(!is.nan(rowSums(x))) else as.integer(!is.nan(x))
rci <- function(x) if(is.matrix(x)) as.integer(rowSums(x) > 0) else as.integer(x > 0)
byaa<-function(x) unlist(tapply(x,ctab$aa,cnorm))

byaas <- function(x) do.call(cbind, sapply(names(acnt), function(y) cnormm(x[,ctab$aa %in% y])))
byrcs <- function(x) sapply(names(acnt), function(y) rcm(x[,ordaa %in% y]))
byrci <- function(x) sapply(names(acnt), function(y) rci(x[,aclist[[y]]]))
