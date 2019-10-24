# getPPIs

## Description
`getPPIs` facilitates the compilation of experimentally identified 
protein-protein interactions (PPIs) from [HitPredict](http://hintdb.hgc.jp/htp/).

#### HitPredict
> Protein-protein interactions from IntAct, BioGRID, HPRD, MINT and DIP are 
> combined, annotated and scored. The reliability score is calculated based 
> on the experimental details of each interaction and the sequence, structure and 
> functional annotations of the interacting proteins.

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
1. __getHitPredict__ - facilitates download of HitPredict data and mapping protein identifiers to stable Entrez gene IDs.
2. __getHomoloGene__ - maps genes to homologs in a species of interest (e.g. map human PPIs to mouse).
3. __getMethods__ - annotates the HitPredict detection methods with more human-readable names.
4. __getPPIs__ - wrapper function that does the work getHitPredict, getHomoloGene, and getPPIs.
5. __buildNetwork__ - builds a PPI graph given some genes of interest (currently only works for mouse).

### Usage Examples:
```R
library(getPPIs)

# 1. Download HitPredict database.
hitpredict <- getHitPredict("HitPredict")

# 2. Map all genes to mouse homologs.
hitpredict <- getHomoloGene(hitpredict, taxid=10090)

# 3. Annotate HitPredict data with method names.
hitpredict <- getMethods(hitpredict)

# 4. Equivalently, steps 1-3 are performed by the wrapper function `getPPIs`:
ppis <- getPPIs("HitPredict", taxid = 10090)

# 5. Given some genes, build a PPI graph:
mygenes <- unique(c(ppis$osEntrezA,ppis$osEntrezB))
g <- buildNetwork(ppis, mygenes, mytaxid=10090) # This is the mouse interactome.

```
See additional examples [here](./examples.R).

## Additional Datasets
`getPPIs` contains several useful datasets, access them in R with the `data()` function:

* `data(musInteractome)` All HitPredict interactions mapped to mouse.
* `data(Wrp)` Wrp-BioID from Spence __et al.__, 2019 <sup>[3]<sup>
* `data(iPSD)` iPSD-BioID from Uezu __et al.__, 2016 <sup>[4]<sup>
* `data(ePSD)` PSD95-BioID from Uezu __et al.__, 2016 <sup>[4]<sup>
* `data(compiled_iPSD)` An iPSD proteome compiled from several studies <sup>[5-n]<sup>

Note: You can download the compiled mouse interactome data
[here](https://github.com/twesleyb/getPPIs/blob/master/data/musInteractome.zip).

## Dependencies
The script `installDBs` automates the installation of several gene ID mapping databases
which are utilized to map organism specific Uniprot entries to Entrez gene
identifiers, a unique stable identifier for every gene in every organism.  

Other libraries utilized by `getPPIs` include:
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
__[1]__ Patil, A. & Nakamura, H. Filtering high-throughput protein-protein interaction data using a combination of genomic features. BMC Bioinformatics 6, 100 (2005).  
__[2]__ Patil, A., Nakai, K. & Nakamura, H. HitPredict: a database of quality assessed protein-protein interactions in nine species. Nucleic Acids Res. 39, D744-9 (2011).  
__[3]__ López, Y., Nakai, K. & Patil, A. HitPredict version 4: comprehensive reliability scoring of physical protein-protein interactions from more than 100 species. Database (Oxford) 2015, (2015).  
__[4]__ Uezu, A. et al. Identification of an elaborate complex mediating postsynaptic inhibition. Science 353, 1123–1129 (2016).  
__[5]__ Spence, E. F. et al. In vivo proximity proteomics of nascent synapses reveals a novel regulator of cytoskeleton-mediated synaptic maturation. Nat. Commun. 10, 386 (2019).  
