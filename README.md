# getPPIs

An R package that facilitates the compilation of experimentally identified, and 
scored interactions (PPIs) from the [HitPredict](http://hintdb.hgc.jp/htp/) database.

#### HitPredict:
> Protein-protein interactions from IntAct, BioGRID, HPRD, MINT and DIP are 
> combined, annotated and scored. The reliability score is calculated based 
> on the experimental details of each interaction and the sequence, structure and 
> functional annotations of the interacting proteins.

For more information about HitPredict, please see authors' original publications:
* Lopez _et al.,_ [Database 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4691340/)
* Patil _et al.,_ [Nucleic Acids Res. 2011](https://www.ncbi.nlm.nih.gov/pubmed/20947562)
* Patil and Nakamura [BMC Bioinformatics 2005](https://www.ncbi.nlm.nih.gov/pubmed/15833142)

## Installation
In R, download the package from GitHub:

```R
devtools::install_github("twesleyb/getPPIs")
```

## Usage
The package contains several key functions:
1. getHitPredict() - facilitates download of HitPredict data and mapping protein identifiers to stable Entrez gene IDs.
2. getHomoloGene() - automates download of NCBI homology database and maps genes to their homolog in a species of interest.
3. getMethods() - annotates the HitPredict data with more human-readable names cooresponding to detection methods.
4. getPPIs() - wrapper function that does the work getHitPredict, getHomoloGene, and getPPIs.
5. buildNetwork() - builds a protein-protein interaction graph given some genes of interest.

To use `getPPIs`, just try this:
```
library(getPPIs)
data(iPSD)
ppis <- getPPIs(organism="HitPredict", iPSD, taxid=10090)

```

Then, build a PPI graph:
```
g <- buildNetwork(ppis,taxid=10090)
```
## Additional Datasets

#### Mouse interactome
PPIs among mouse proteins were compiled. You can easily load this data with:
```
data(musInteractome)

```

#### BioID datasets
The package contains several published iBioID datasets. To access them, use the
`data()` function:

```
library(getPPIs)

data(Wrp)
data(iPSD)
data(compiled_iPSD)
data(ePSD)
```
