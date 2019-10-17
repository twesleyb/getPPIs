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
data(iPSD)

# Get PPIs amongst iPSD proteins (mouse).
# Keep all PPIs that are homologous to mouse genes.
hitpredict <- getPPIs(organism="HitPredict", mygenes=iPSD, taxid=10090)

# Build a protein-protein interaction network.


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


