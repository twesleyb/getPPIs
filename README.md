# getPPIs
Compiling experimentally identified protein-protein interactions (PPIs) from 
the __HitPredict__ database.

## Motivation
This package facilitates downloading and compiling PPIs from the HitPredict database.

## HitPredict
[HitPredict](http://hintdb.hgc.jp/htp/)
> Protein-protein interactions from IntAct, BioGRID, HPRD, MINT and DIP are combined, annotated and scored. The reliability score is calculated based on the experimental details of each interaction and the sequence, structure and functional annotations of the interacting proteins

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

