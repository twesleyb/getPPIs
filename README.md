# getPPIs
Compiling experimentally identified protein-protein interactions (PPIs) from 
the __HitPredict__ database.

## Motivation
This package facilitates downloading and compiling PPIs from the HitPredict database.

## HitPredict
[HitPredict](http://hintdb.hgc.jp/htp/)
> Protein-protein interactions from IntAct, BioGRID, HPRD, MINT and DIP are combined, annotated and scored. The reliability score is calculated based on the experimental details of each interaction and the sequence, structure and functional annotations of the interacting proteins

For more information, please see authors' original publications:
* Lopez _et al.,_ [Database 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4691340/)
* Patil _et al.,_ [Nucleic Acids Res. 2011](https://www.ncbi.nlm.nih.gov/pubmed/20947562)
* Patil and Nakamura [BMC Bioinformatics 2005](https://www.ncbi.nlm.nih.gov/pubmed/15833142)

## Usage
Download the getPPIs R package from GitHub.
```
devtools::install_github("twesleyb/getPPIs")

library(getPPIs)

```

