#' KEGG Orthology (KO) annotations for ribosomal genes.
#'
#' @format A character vector.
#'
RPKOs <- c(
    "K02945","K02967","K02982","K02986","K02988","K02990","K02992","K02994",
    "K02996","K02946","K02948","K02950","K02952","K02954","K02956","K02959",
    "K02961","K02963","K02965","K02968","K02970","K02981","K02985","K02984",
    "K02987","K02989","K02991","K02993","K02995","K02997","K02947","K02949",
    "K02951","K02953","K02955","K02958","K02957","K02960","K02962","K02964",
    "K02966","K02969","K02971","K02973","K02974","K02975","K02976","K02978",
    "K02977","K02979","K02980","K02983","K02998","K02863","K02886","K02906",
    "K02926","K02931","K02933","K02935","K02939","K02864","K02867","K02869",
    "K02871","K02874","K02876","K02878","K02879","K02881","K02884","K02887",
    "K02888","K02890","K02892","K02895","K02897","K02899","K02902","K02904",
    "K02907","K02909","K02911","K02913","K02914","K02916","K02919","K07590",
    "K02925","K02930","K02932","K02934","K02937","K02936","K02938","K02940",
    "K02866","K02865","K02868","K02870","K02873","K02872","K02875","K02877",
    "K02880","K02883","K02882","K02885","K02889","K02891","K02894","K02893",
    "K02896","K02898","K02901","K02900","K02903","K02905","K02908","K02910",
    "K02912","K02915","K02918","K02917","K02920","K02922","K02921","K02923",
    "K02924","K02927","K02928","K02929","K02941","K02942","K02943","K02944",
    "K01977","K01980","K01985","K01979","K01982","K01981"
)

#' Codon usage in healthy human gut microbiome.
#'
#' A \code{codonTable} object with codon counts for sequences of the human gut
#' metagenome, from a healthy individual. Raw sequences are from
#' \href{https://go.nature.com/2G7QdpC}{Quin et al. 2014},
#' processed, assembled and annotated (KEGG Orthology) as descrbed in
#' \href{https://bit.ly/2DRfiz6}{Fabijanic and Vlahovicek 2016}.
#'
#' @format A \code{codonTable} object.
#'
#' @source
#' \href{https://go.nature.com/2G7QdpC}{Quin et al. 2014};
#' \href{https://bit.ly/2DRfiz6}{Fabijanic and Vlahovicek 2016}
"HD59"

#' Codon usage in human gut microbiome in liver cirrhosis.
#'
#' A \code{codonTable} object with codon counts for sequences of the human gut
#' metagenome, from an individual with liver cirrhosis. Raw sequences are from
#' \href{https://go.nature.com/2G7QdpC}{Quin et al. 2014},
#' processed, assembled and annotated (KEGG Orthology) as descrbed in
#' \href{https://bit.ly/2DRfiz6}{Fabijanic and Vlahovicek 2016}.
#'
#' @format A \code{codonTable} object.
#'
#' @source
#' \href{https://go.nature.com/2G7QdpC}{Quin et al. 2014};
#' \href{https://bit.ly/2DRfiz6}{Fabijanic and Vlahovicek 2016}
"LD94"

#' Codon usage based KO enrichment analysis results
#' from a healthy human gut microbiome.
#'
#' @format An \code{AnnotatedDataFrame} object.
#' See `?enrichment` for description.
#'
#' @source
#' \href{https://go.nature.com/2G7QdpC}{Quin et al. 2014};
#' \href{https://bit.ly/2DRfiz6}{Fabijanic and Vlahovicek 2016}
"HD59_KO"


#' Codon usage based KO enrichment analysis results
#' from an gut microbiome of an individual with liver cirrhosis.
#'
#' @format An \code{AnnotatedDataFrame} object.
#' See `?enrichment` for description.
#'
#' @source
#' \href{https://go.nature.com/2G7QdpC}{Quin et al. 2014};
#' \href{https://bit.ly/2DRfiz6}{Fabijanic and Vlahovicek 2016}
"LD94_KO"

#' Codon usage based KEGG Pathway enrichment analysis results
#' from a healthy human gut microbiome.
#'
#' @format An \code{AnnotatedDataFrame} object.
#' See `?enrichment` for description.
#'
#' @source
#' \href{https://go.nature.com/2G7QdpC}{Quin et al. 2014};
#' \href{https://bit.ly/2DRfiz6}{Fabijanic and Vlahovicek 2016}
"HD59_PATHWAYS"


#' Codon usage based  KEGG Pathway enrichment analysis results
#' from an gut microbiome of an individual with liver cirrhosis.
#'
#' @format An \code{AnnotatedDataFrame} object.
#' See `?enrichment` for description.
#'
#' @source
#' \href{https://www.nature.com/articles/nature13568}{Quin et al. 2014};
#' \href{https://bit.ly/2DRfiz6}{Fabijanic and Vlahovicek 2016}
"LD94_PATHWAYS"
