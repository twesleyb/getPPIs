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
See additional usage examples [here](./examples.R).

## Additional Datasets
`getPPIs` contains several useful datasets, access them in R with the `data()` function:

* `data(musInteractome)` All HitPredict interactions mapped to mouse.
* `data(Wrp)` Wrp-BioID from Spence _et al._, 2019.<sup>[3]<sup>
* `data(iPSD)` iPSD-BioID from Uezu _et al._, 2016.<sup>[4]<sup>
* `data(ePSD)` PSD95-BioID from Uezu _et al._, 2016.<sup>[4]<sup>
* `data(compiled_iPSD)` An iPSD proteome compiled from several studies.<sup>[4-8]<sup>

Download the compiled mouse interactome
[here](https://github.com/twesleyb/getPPIs/blob/master/data/musInteractome.zip).

## The Compiled iPSD proteome
The proteome of the inhibitory post-synaptic density (iPSD) is understudied relative to its excitatory, PSD conterpart. 
The compiled iPSD is a list of proteins that have been identified at inhibitory synapses by a variety of methods.

| Study                   | Method                                                                         |  nProteins  | Ref.|
|:-----------------------:|:------------------------------------------------------------------------------:|:-----------:|:---:|
| Heller _et al._, 2012   | Affinity purification of Venus eGFP-tagged GABAARa1 transgenic mouse           | 18  | [5] |
| Kang _et al._, 2014     | Affinity purification of His6-FLAG-YFP-NL2 transgenic mouse                    | 75  | [6] |
| Uezu _et al._, 2016     | Gephyrin, InSyn1, and Collybistin iBioID                                       | 181 | [4] |
| Loh _et al._, 2016      | Lenti-viral expression of Slitrk3 and Nlgn2 APEX in cultured rat neurons       | 42  | [7] |
| Nakamura _et al._, 2016 | Affinity purification of pHluorin tagged GABAARa2 transgenic mouse             | 174 | [8] |
| Yamasaki _et al._, 2017 | GABAARa1 immunoprecipitation                                                   | 7   | [7] |

You can download the compiled iPSD proteome [here](./link).

## Dependencies
The script `installDBs` automates the installation of several gene ID mapping databases
which are utilized to map organism specific Uniprot entries to Entrez gene IDs.

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

[1] Patil, A. & Nakamura, H. Filtering high-throughput protein-protein interaction data using a combination of genomic features. BMC Bioinformatics 6, 100 (2005).  

[2] Patil, A., Nakai, K. & Nakamura, H. HitPredict: a database of quality assessed protein-protein interactions in nine species. Nucleic Acids Res. 39, D744-9 (2011).  

[3] López, Y., Nakai, K. & Patil, A. HitPredict version 4: comprehensive reliability scoring of physical protein-protein interactions from more than 100 species. Database (Oxford) 2015, (2015).  

[4] Uezu, A. et al. Identification of an elaborate complex mediating postsynaptic inhibition. Science 353, 1123–1129 (2016).  

[5] Spence, E. F. et al. In vivo proximity proteomics of nascent synapses reveals a novel regulator of cytoskeleton-mediated synaptic maturation. Nat. Commun. 10, 386 (2019).  

## License
This program is free software; you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. 
If not, see http://www.gnu.org/licenses/.
