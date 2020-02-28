# getPPIs

## Description
_getPPIs_ is an __R__ package that facilitates the compilation of protein-protein 
interactions from [HitPredict](http://hintdb.hgc.jp/htp/).

#### HitPredict
> HitPredict is a database of experimentally identified, confidence scored 
> protein-protein interactions (PPIs) compiled from several databases, including:
> [IntAct](https://www.ebi.ac.uk/intact/), 
> [BioGRID](https://thebiogrid.org/), 
> [HPRD](https://hprd.org/), 
> [MINT](https://mint.bio.uniroma2.it/), and 
> [DIP](https://dip.doe-mbi.ucla.edu/dip/Main.cgi).

For more information about HitPredict, please see authors' original publications:
* [Lopez _et al.,_ 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4691340/) <sup>[1]<sup>
* [Patil _et al.,_ 2011](https://www.ncbi.nlm.nih.gov/pubmed/20947562) <sup>[2]<sup>
* [Patil and Nakamura 2005](https://www.ncbi.nlm.nih.gov/pubmed/15833142) <sup>[3]<sup>

## Installation
In R, download `getPPIs` from GitHub:

```R
devtools::install_github("twesleyb/getPPIs")
```

## Usage
The package contains several key functions:
1. __getPPIs__ - Build a protein-protein interaction database. 
A wrapper function that the work getHitPredict, getHomoloGene, and getInteractionMethods.
2. __mapIDs__ - Map gene IDs from one format to another (e.g. uniprot to entrez).
3. __getHomologs__ - Map gene identifiers to homologs in another species.
4. __buildNetwork__ - builds a PPI graph given some genes of interest (currently only works for mouse).

### Usage Examples:
```R
library(getPPIs)

# 1. Download HitPredict database.
hitpredict <- getHitPredict("HitPredict")

# 2. Map all genes to mouse homologs.
hitpredict <- getHomoloGene(hitpredict, species="mouse")

# 3. Annotate HitPredict data with method names.
hitpredict <- getInteractionMethods(hitpredict)

# 4. Equivalently, steps 1-3 are performed by the wrapper function `getPPIs`:
ppis <- getPPIs("HitPredict", species = "mouse")

# 5. Given some genes, build a PPI graph:
mygenes <- unique(c(ppis$osEntrezA,ppis$osEntrezB))
g <- buildNetwork(ppis, mygenes, mytaxid=10090) # This is the mouse interactome.

```
See additional usage examples [here](./examples.R).

## Dependencies
The script `installDBs` automates the installation of several gene ID mapping databases
which are utilized to map organism specific Uniprot entries to Entrez gene IDs.

Other key libraries utilized by `getPPIs` include:
* data.table
* dplyr
* igraph
* xml2
* rvest
* tools
* AnnotationDbi
* ontologyIndex
* BiocManager

## References

[1] Patil, A. & Nakamura, H. Filtering high-throughput protein-protein interaction data using a combination of genomic features. BMC Bioinformatics 6, 100 (2005).  

[2] Patil, A., Nakai, K. & Nakamura, H. HitPredict: a database of quality assessed protein-protein interactions in nine species. Nucleic Acids Res. 39, D744-9 (2011).  

[3] López, Y., Nakai, K. & Patil, A. HitPredict version 4: comprehensive reliability scoring of physical protein-protein interactions from more than 100 species. Database (Oxford) 2015, (2015).  

## License
This program is free software; you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. 
If not, see http://www.gnu.org/licenses/.
