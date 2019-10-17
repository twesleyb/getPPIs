#------------------------------------------------------------------------------
## Installation.
#------------------------------------------------------------------------------

# Install from github.
devtools::install_github("twesleyb/getPPIs")

#------------------------------------------------------------------------------
## Building an interactome from scratch.
#------------------------------------------------------------------------------

# Load the library.
library(getPPIs)

# Download all PPIs in the HitPredict database.
# This function will also annotate Genes with Entrez IDs, and discard those 
# which are not mapped to Entrez..
hitpredict <- getHitPredict(organism="HitPredict")

# Map genes to mouse homologs.
hitpredict <- getHomoloGene(hitpredict, taxid = 10090)

# Annotate hitpredict data with more descriptive method names.
# Keep all experimental evidence, no confidence score cutoff.
hitpredict <- getMethods(hitpredict)

#------------------------------------------------------------------------------
## Using the pre-built mouse Interactome.
#------------------------------------------------------------------------------

library(getPPIs)

# Load the mouse interactome.
data(musInteractome)
data(compiled_iPSD)

# Build Synaptosome PPI graph.
g <- buildNetwork(hitpredict=musInteractome, mygenes=compiled_iPSD, taxid=10090)
