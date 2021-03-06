% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codonUsage-visualization.R
\docType{methods}
\name{intraBplot}
\alias{intraBplot}
\alias{intraBplot,codonTable,codonTable-method}
\title{Intra-samples Karlin B plot}
\usage{
intraBplot(x, y, names = c("x", "y"), variable, ribosomal = FALSE,
  size = 1, alpha = 0.5)

\S4method{intraBplot}{codonTable,codonTable}(x, y, names = c("x", "y"),
  variable, ribosomal = FALSE, size = 1, alpha = 0.5)
}
\arguments{
\item{x, y}{Objects of \code{codonTable} class.}

\item{names}{Character vector of length 2, giving names for samples.}

\item{variable}{A character, name of the function that will be used
to calculate CU statistic values for plotting. Must be one of
the following: \code{c("MILC", "B", "MCB", "ENCprime")}.}

\item{ribosomal}{Logical, whether to indicate ribosomal genes in the plot.
Default is \code{FALSE}.}

\item{size}{Numeric, indicating points' size}

\item{alpha}{Numeric, between 0 and 1, indicating points' transparency
(default is 0.1).}
}
\value{
A \code{ggplot} object.
}
\description{
Plot CU frequency distances between two samples
(given by \code{x} and \code{y}).
}
\examples{
require(ggplot2)
# calculate MILC distance to the average CU of the example DNA sequences,
# and to the average CU of ribosomal genes among the example DNA sequences
milc <- MILC(LD94, self = TRUE, ribosomal = TRUE)

intraBplot(x = HD59, y = LD94, names = c("HD59", "LD94"),
           variable = "MILC", size = 3)

}
