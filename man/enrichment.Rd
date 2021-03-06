% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment.R
\docType{methods}
\name{enrichment}
\alias{enrichment}
\alias{enrichment,crossTab-method}
\title{Enrichment analysis for codon usage (CU) data.}
\usage{
enrichment(x, pvalueCutoff = numeric(), pAdjustMethod = "BH",
  padjCutoff = numeric())

\S4method{enrichment}{crossTab}(x, pvalueCutoff = numeric(),
  pAdjustMethod = "BH", padjCutoff = numeric())
}
\arguments{
\item{x}{A \code{crossTab} object}

\item{pvalueCutoff}{Numeric, discard categories with p value below this
threshold. By default, no threshold is set (\code{numeric()}).}

\item{pAdjustMethod}{Character, one of the \code{p.adjust.methods}.}

\item{padjCutoff}{Numeric, discard categories with adjusted p value below
this threshold. By default, no threshold is set (\code{numeric()}).}
}
\value{
An \code{AnnotatedDataFrame} object, or a list of those; data in
each object has category values in rows, and the following columns:
   \itemize{
     \item category, a character vector of annotation categories
     \item all, a numeric vector of integers, coresponding to sequence
         counts for each annotation category, in the background gene set
         (universe).
     \item a numeric vector(s) of integers, coresponding to sequence counts
         for each annotation category, in the set of genes for which
         enrichment is calculated, i.e. the predefined subset of (usually
         highly expressed) genes in the universe (named for the
         corresponding `crossTab` column).
     \item enrichment, calculated as the ratio: (scaled sample counts -
         scaled backg. counts) / scaled backg. counts * 100,
         where scaling means that sample counts are simply increased by 1,
         and background counts are multiplied by ratio of summed sample
         counts and summed backgroun counts, and also increased by 1
     \item M, log ratios of scaled counts
     \item A, mean average of scaled counts
     \item pvals, p values for exact binomial test
     \item padj, p values corrected by BH method.
   }
}
\description{
Performs enrichment analysis, given a contongency table of codon counts.
p values are calculated by binomial test, adjustment for multiple testing
can be performed by any of the \code{p.adjust.methods}.
}
\examples{
require(Biobase)

# create contingency table
s <- getKO(HD59)
v <- as.numeric(MELP(HD59, ribosomal = TRUE))
ct <- crossTab(s, v)

# enrichment analysis
enr <- enrichment(ct)
enr # for help, see `?Biobase::AnnotatedDataFrame`
head(pData(enr))

enr <- enrichment(ct, pAdjustMethod = "holm")
head(pData(enr))

enr <- enrichment(ct, pvalueCutoff = 0.05)
head(pData(enr))

enr <- enrichment(ct, padjCutoff = 0.05)
head(pData(enr))

}
