# getPPIs

## Description
_getPPIs_ is an __R__ package containing a compilation of experimentally identified,
confidence scored protein-protein interactions (PPIs) from
[HitPredict](http://hintdb.hgc.jp/htp/).  Simply load a `data.frame` containing the PPI data using
`data(musInteractome)` or `data(hsInteractome)`.

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

```R
library(getPPIs)

# load the HitPredict dataset
data(hsInteractome)

# load the HitPredict dataset, mapped to mouse genes
data(musInteractome)

```


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
