% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codonUsage-visualization.R
\docType{methods}
\name{Bplot}
\alias{Bplot}
\alias{Bplot,character,character,matrix-method}
\alias{Bplot,numeric,numeric,missing-method}
\title{Karlin B plot}
\usage{
Bplot(x, y, data, annotations = character(), ribosomal = FALSE,
  reference = list(), size = 1, alpha = 0.5)

\S4method{Bplot}{character,character,matrix}(x, y, data,
  annotations = character(), ribosomal = FALSE, reference = list(),
  size = 1, alpha = 0.5)

\S4method{Bplot}{numeric,numeric,missing}(x, y, data,
  annotations = character(), ribosomal = FALSE, reference = list(),
  size = 1, alpha = 0.5)
}
\arguments{
\item{x, y}{Character, both must be in \code{colnames(data)}, or numeric
vectors of CU statistic values for two subsets of genes. If numeric,
the vectors must be of the same length.}

\item{data}{A matrix with CU statistic values for subsets of genes
in columns.}

\item{annotations}{A character vector giving KO annotations for sequences
for which the CU values were calculated, must be of length
\code{nrow(data)}.}

\item{ribosomal}{Logical, whether to indicate ribosomal genes in the plot.
Default is \code{FALSE}, if set to \code{TRUE}, then \code{annotation}
must be given.}

\item{reference}{A named list of length 1, containing either a logical
vector of \code{nrow(data)} of reference genes to be indicated
on the plot, or a character vector (of any length) of the reference
genes' anotations. If latter is the case, then \code{annotation}
must be given.}

\item{size}{Numeric, indicating points' size}

\item{alpha}{Numeric, between 0 and 1, indicating points' transparency
(default is 0.1).}
}
\value{
A \code{ggplot} object.
}
\description{
Plot distances of each gene's CU frequency to specified gene (sub)sets
(given by \code{x} and \code{y}).
}
\examples{
require(ggplot2)

# calculate MILC distance to the average CU of the example DNA sequences,
# and to the average CU of ribosomal genes among the example DNA sequences
milc <- MILC(LD94, self = TRUE, ribosomal = TRUE)

Bplot(x = "ribosomal", y = "self", data = milc,
      ribosomal = TRUE, annotations = getKO(LD94),
      size = 3) +
    labs(x = "MILC distance to ribosomal genes",
         y = "MILC distance to genes' average CU")

}
