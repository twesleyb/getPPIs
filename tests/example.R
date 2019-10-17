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
ppis <- getPPIs(organism="HitPredict", mygenes=iPSD, taxid=10090)

# Build a ppi network.
library(org.Mm.eg.db)

sif <- ppis %>% dplyr::select(osEntrezA, osEntrezB)

nodes <- unique(c(sif$osEntrezA, sif$osEntrezB))

symbols <- AnnotationDbi::mapIds(org.Mm.eg.db,
				 keys = as.character(nodes),
				 column = "SYMBOL",
				 keytype = "ENTREZID",
				 multiVals = "first")
noa <- data.table::data.table(node = nodes, symbol = symbols)

# Check.
if (sum(is.na(symbols)) == 0) {
	message("All node Entrez IDs mapped to gene symbols!")
}

# Build igraph object.
library(igraph)
g <- graph_from_data_frame(sif, directed = FALSE, vertices = nodes)

# Insure that the graph is simple.
g <- simplify(g)

# You can save these files to work with them in Cystoscape.
data.table::fwrite(noa,"noa.csv")
data.table::fwrite(sif,"sif.csv")

# Alternatively, check out the RCy3 package to interface with Cytscape from R!

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


