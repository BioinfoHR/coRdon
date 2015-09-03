make.contable <- function(melpSet, variable = "melp", category = "KO", thresholds = c(1), percentiles = NULL) {
  
  all <- as.vector(table(melpSet[[category]]))
  result <- data.frame(all = all)
  
  if(!is.null(percentiles)) {
    top.perc <- sapply(percentiles, function(x) {
      as.vector(table(melpSet[melpSet[[variable]] >= quantile(melpSet[[variable]], x), category]))
    })
    top.perc <- as.data.frame(top.perc)
    names(top.perc) <- make.names(paste("top", perc, sep="_"))
    
    result <- cbind(result, top.perc)
  } else {
    top.perc <- NULL
  }
  
  if(!is.null(thresholds)) {
    top.thresh <- sapply(thresholds, function(x) {
      as.vector(table(melpSet[melpSet[[variable]] >= x, category]))
    })
    top.thresh <- as.data.frame(top.thresh)
    names(top.thresh) <- make.names(paste("gt", thresholds, sep="_"))
    
    result <- cbind(result, top.thresh)
  } else {
    top.thresh <- NULL
  }
  
  categories <- levels(melpSet[[category]])
  data.frame(result, category = categories, row.names = categories)
}