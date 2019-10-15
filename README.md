# getPPIs
An R package that facilitates the compilation of experimentally identified
protein-protein interactions (PPIs) from the [HitPredict](http://hintdb.hgc.jp/htp/)
database.

## Installation
Download the package from GitHub with devtools.

```
devtools::install_github("twesleyb/getPPIs")
library(getPPIs)
```

## HitPredict
> Protein-protein interactions from IntAct, BioGRID, HPRD, MINT and DIP are combined, annotated and scored. The reliability score is calculated based on the experimental details of each interaction and the sequence, structure and functional annotations of the interacting proteins

For more information, please see authors' original publications:
* Lopez _et al.,_ [Database 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4691340/)
* Patil _et al.,_ [Nucleic Acids Res. 2011](https://www.ncbi.nlm.nih.gov/pubmed/20947562)
* Patil and Nakamura [BMC Bioinformatics 2005](https://www.ncbi.nlm.nih.gov/pubmed/15833142)

## Approach

### Building an interactome from scratch:
#### getPPIs()
The function `getPPIs()` automates the compilation of experimentally determined PPIs. It is a wrapper around the following functions:
1. __getHitPredict()__: HitPredict database is downloaded and genes are mapped to stable identifiers (Entrez IDS).
2. __getHomoloGene()__: Genes are mapped to homologs in the NCBI database.
3. __getMethods()__: Detection methods are annotated.


## Mouse interactome

```
Mapping Uniprot IDs to Entrez IDs...
   ... 94.45 % of yeast Uniprot IDs were successfully mapped to Entrez IDs.
   ... 99.27 % of worm Uniprot IDs were successfully mapped to Entrez IDs.
   ... 97.57 % of fly Uniprot IDs were successfully mapped to Entrez IDs.
   ... 51.43 % of zebrafish Uniprot IDs were successfully mapped to Entrez IDs.
   ... 97.4 % of frog Uniprot IDs were successfully mapped to Entrez IDs.
   ... 93.31 % of chicken Uniprot IDs were successfully mapped to Entrez IDs.
   ... 98.47 % of human Uniprot IDs were successfully mapped to Entrez IDs.
   ... 67.95 % of dog Uniprot IDs were successfully mapped to Entrez IDs.
   ... 89.29 % of pig Uniprot IDs were successfully mapped to Entrez IDs.
   ... 88.29 % of bovine Uniprot IDs were successfully mapped to Entrez IDs.
   ... 96.7 % of mouse Uniprot IDs were successfully mapped to Entrez IDs.
   ... 95.73 % of rat Uniprot IDs were successfully mapped to Entrez IDs.
17,761 (3.54%) of all interactions (rows) are incompletely mapped to Entrez IDs and will be discarded.
483,477 Protein-protein interactions from 12 species were compiled from the HitPredict database and successfully mapped to Entrez IDs.
```

#### Summary
21.95% of HitPredict interactions were mapped to homologous genes in mouse.  
377,352 protein-protein interactions among 15,945 mouse genes were compiled from 10 species.  


## BioID datasets
The package contains several published iBioID datasets. To access them:

```
data(wrp)
data(iPSD)
data(compiled_iPSD)
data(ePSD)
```

