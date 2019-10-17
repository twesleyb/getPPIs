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
## Using the pre-built mouse Interactome.
#------------------------------------------------------------------------------

library(getPPIs)

# Load the mouse interactome.
data(musInteractome)
data(compiled_iPSD)

# Build Synaptosome PPI graph.
# FIXME: buildNetwork to supress unwanted text output!
g <- buildNetwork(hitpredict=musInteractome, mygenes=compiled_iPSD, taxid=10090)

# Community detection with the Leiden algorithm.
library(leiden)

partition <- leiden(g, partition_type = "ModularityVertexPartition",
		    resolution_parameter = 1, seed = NULL, n_iterations = -1)

nodes <- names(V(g))
names(nodes) <- vertex_attr(g, "symbol")

"insyn2" %in% tolower(names(nodes))

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


