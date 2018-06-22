##coRdon
### Codon Usage Analysis and Prediction of Gene Expressivity
<br>

R package for analysis and visualization of codon usage in DNA sequences.

Main functionalities:
* calculates different measures of CU bias and CU-based predictors 
of gene expressivity
* performs gene set enrichment analysis for unannotated or KEGG/COG 
annotated DNA sequences
* implements several methods for visualization of codon usage 
and enrichment analysis results.

Although specifically aimed at the analysis of metagenomic samples, 
coRdon allows for the inspection and quantification of codon usage 
in DNA sequences using any of the 20 different variants of genetic code,
with the additional options to include stop codons and alternative start 
codons in calculations.  

The following statistics are implemented in the package:

* ENC, effective number of codons 
([Wright, 1990](https://www.ncbi.nlm.nih.gov/pubmed/2110097)),  
* its modified version ENC' 
([Novembre, 2002](https://www.ncbi.nlm.nih.gov/pubmed/12140252)),  
* a measure of codon bias, termed B 
([Karlin and Mrazek, 1996](https://www.ncbi.nlm.nih.gov/pubmed/11489855)),  
* and related measure of expression, E
([Karlin and Mrazek, 2000](https://www.ncbi.nlm.nih.gov/pubmed/10960111)),  
* maximum likelihood codon bias, MCB 
([Urrutia and Hurst, 2001](https://www.ncbi.nlm.nih.gov/pubmed/2110097)),  
* MILC, Measure Independent of Length and Composition, and  
* MELP, MILC-based Expression Level Predictor
([Supek and Vlahovicek, 2005](https://www.ncbi.nlm.nih.gov/pubmed/16029499)),
* SCUO, synonymous codon usage orderliness 
([Wan et al., 2004](https://www.ncbi.nlm.nih.gov/pubmed/15222899)),
* Codon Adaptation Index, CAI
([Sharp and Li, 1987](https://www.ncbi.nlm.nih.gov/pubmed/3547335)),  
* frequency of optimal codons, Fop
([Ikemura 1981](https://www.ncbi.nlm.nih.gov/pubmed/6175758)),  
* gene codon bias, GCB
([Merkl, 2003](http://www.ncbi.nlm.nih.gov/pubmed/14708578)).

The package also implements B plot for visualization of CU bias, 
both within a single sample and between different samples for which 
CU bias statistics are calculated.  

Additionally, if the input sequences are annotated in either KEG or COG 
orthology database, functional analysis can be performed in order to determine
significantly enriched functions in the imput sample. This is aimed 
particularly at metagenomic samples, as a way of determining 
functional fingerprint of a microbial community.
There are also several methods for visualisation of enrichment analysis results,
including MA-like plot and bar plot.  

***

## Geting started

To install coRdon, run the following in R:
```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("coRdon")
```

The developmental version can be installed directly from from GitHub: 
```{r}
devtools::install_github("BioinfoHR/coRdon")
```

For worked example on how to do analysis of codon usage with coRdon, 
please see the packaage vignette.

***

<i><img src="http://bioinfo.hr/wp-content/themes/theme1414/images/logo.png" alt="bioinfo.hr" title=""></i>