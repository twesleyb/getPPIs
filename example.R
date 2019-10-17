#------------------------------------------------------------------------------
## Installation.
#------------------------------------------------------------------------------

# Install from github.
devtools::install_github("twesleyb/getPPIs")

#------------------------------------------------------------------------------
## Download, compile, annotate, and subset PPIs in a single step.
#------------------------------------------------------------------------------

# Load the library.
library(getPPIs)

# Load an example dataset.
data(compiled_iPSD)

# Get all HitPredict PPIs. 
# Keep all PPIs that are homologous to mouse genes.
ppis <- getPPIs(organism="HitPredict", taxid=10090)

#------------------------------------------------------------------------------
## Using the prebuild mouse Interactome.
#------------------------------------------------------------------------------

library(getPPIs)

# Load the mouse interactome.
data(musInteractome)
data(compiled_iPSD)

# Build Synaptosome PPI graph.
g <- buildNetwork(hitpredict=musInteractome, mygenes=compiled_iPSD, taxid=10090)

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


