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
# This function will also annotate HitPredict proteins with Entrez IDs, a 
# stable, unique gene identifier, and discard those which are not mapped.
hitpredict <- getHitPredict(organism="HitPredict")

# Map genes to homologs in your species of interest (mouse = 10090)..
hitpredict <- getHomoloGene(hitpredict, taxid = 10090)

# Annotate hitpredict data with more descriptive method names.
# Keep all experimental evidence, no confidence score cutoff.
hitpredict <- getMethods(hitpredict)

# getPPIs does the previous three steps.
hitpredict <- getPPIs("HitPredict", taxid=10090)

#------------------------------------------------------------------------------
## Using the pre-built mouse Interactome.
#------------------------------------------------------------------------------
# If you work in mice, and you just want to find PPIs among your proteins of
# interest, you can skip building your own interactome, and just load the 
# pre-built mouse interactome.

library(getPPIs)

# Load the mouse interactome.
data(musInteractome)

# Load your genes of interest into R.
# For example, Uniprot identifiers for the iPSD proteome were mapped to Entrez
# IDs using the MGI website: http://www.informatics.jax.org/batch.
# Your data just needs to have a column containing the word 'entrez' (case-insensitive).
# Data file should be `.txt` or `.csv` format. 
myfile <- file.path("./data","MGIBatchReport_20191024_171834.txt")

# Collect PPIs among your proteins and build a PPI graph.
ppis <- buildNetwork(musInteractome,mygenes=myfile,taxid=10090)

#------------------------------------------------------------------------------
## Using the BioID datasets.
#------------------------------------------------------------------------------

library(getPPIs)

# Load the data.
data(iPSD)

data(Wrp)
data(ePSD)

# Build some graphs.
ipsd <- getPPIs::buildNetwork(hitpredict=musInteractome, mygenes=iPSD, taxid=10090)

wrp <- buildNetwork(hitpredict=musInteractome, mygenes=Wrp, taxid=10090)
epsd <- buildNetwork(hitpredict=musInteractome, mygenes=ePSD, taxid=10090)

#------------------------------------------------------------------------------
## Mapping gene identifiers programmatically.
#------------------------------------------------------------------------------
# Load the InSyn1 GFP-trap data from Uezu et al., 2016. Map gene names to 
# Entrez Ids.

library(getPPIs)
library(readxl)

## Load your data with readxl::read_excel.
# File path to your data.
myfile <- file.path("./data","Uezu_et_al_2016_TableS6.xlsx")
data <- read_excel(myfile)
genes <- data$"Gene name"
head(genes)

mapIDs(identifiers=genes,from="symbol",to="entrez",species="mouse")

mapIDs <- function(identifiers,from,to,species,...) {
	# Wrapper around AnnotationDbi::mapIds()
	#c("ACCNUM", "ALIAS", "ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS",
	#  "ENTREZID", "ENZYME", "EVIDENCE", "EVIDENCEALL", "GENENAME", 
	#  "GO", "GOALL", "IPI", "MGI", "ONTOLOGY", "ONTOLOGYALL", "PATH",
	#  "PFAM", "PMID", "PROSITE", "REFSEQ", "SYMBOL", "UNIGENE", "UNIPROT")
	# Imports.
	require(getPPIs)
	# Get organism specific mapping database.
	annotationDBs <- mappingDBs()
	orgDB <- unlist(annotationDBs[sapply(annotationDBs, "[", 3) == tolower(species)])
	names(orgDB) <- sapply(strsplit(names(orgDB),"\\."),"[",2)
	suppressPackageStartupMessages({
		eval(parse(text = paste0("require(", orgDB[["database"]], ",quietly=TRUE)")))
	})
	osDB <- eval(parse(text = orgDB[["database"]]))
	# Map gene identifiers. Suppress messages.
      	suppressMessages({
		output <- AnnotationDbi::mapIds(osDB,
						keys = as.character(identifiers),
						column = toupper(to),
						keytype = keytype,
						multiVals = "first",
						)
	})
	# Check that all nodes (entrez) are mapped to gene symbols.
	not_mapped <- entrez[is.na(output)]
	if (sum(is.na(output)) != 0) {
	       message(paste("Unable to map", length(not_mapped), "Entrez IDs to gene symbols!"))
	}
	names(output) <- indentifiers
	return(output)
}


#------------------------------------------------------------------------------
## Using the RCy3 module to interact with Cytoscape. 
#------------------------------------------------------------------------------

# Install getPPIs
devtools::install_github("twesleyb/getPPIs")

library(getPPIs)

# Load the mouse interactome.
data(musInteractome)
data(compiled_iPSD)

# Build Synaptosome PPI graph.
g <- buildNetwork(musInteractome, compiled_iPSD)

# Node attributes.
nodes <- names(V(g))
names(nodes) <- vertex_attr(g, "symbol")
noa <- data.frame(node=nodes,symbol=names(nodes))
rownames(noa) <- noa$node

# Send to Cytoscape.
library(RCy3)
cytoscapePing()

# Send graph and node attributes to cytoscape.
createNetworkFromIgraph(g,"compiled_iPSD")
loadTableData(noa)

# Apply FD layout.
layoutNetwork(layout.name = "force-directed") 

# Customize the style
library(TBmiscr)
style <- "SynProt"
defaults <- list(NODE_FILL = col2hex("grey"),
                 NODE_SHAPE="Ellipse",
                 NODE_SIZE=55,
                 EDGE_WIDTH=2.0,
                 EDGE_TRANSPARENCY=120)

nodeLabels <- mapVisualProperty('node label','symbol','p')
createVisualStyle(style, defaults, list(nodeLabels))
lockNodeDimensions(TRUE, style)

#------------------------------------------------------------------------------
## Community detection with the Leiden algorithm.
#------------------------------------------------------------------------------

library(leiden)

# Loop to partition network at several resolutions.
profile <- list()
resolution <- seq(0,1,0.01)
for (i in 1:length(resolution)){
	print(i)
	partition <- leiden(g, partition_type = "ModularityVertexPartition",
			    resolution_parameter = resolution[i], seed = NULL, n_iterations = -1)
	profile[[i]] <- partition
}

# Number of clusters
k <- unlist(lapply(profile,function(x) length(unique(x))))

# Loop to calculate modularity.
q <- unlist(lapply(profile, function(x) modularity(g, membership=x)))

q[q==max(q)]


nodes <- names(V(g))
names(nodes) <- vertex_attr(g, "symbol")

# Add node attributes.
g <- set_vertex_attr(g, "p1", value=partition)

vertex_attr_names(g)

#------------------------------------------------------------------------------
## Analysis of the mouse interactome.
#------------------------------------------------------------------------------

# Load the data.
library(getPPIs)
data(musInteractome)

# Build a network.
mygenes <- unique(c(musInteractome$osEntrezA,musInteractome$osEntrezB))

# Add ability to use cuttoff score... 
data <- buildNetwork(musInteractome,mygenes)

g <- data$network
ppis <- data$ppis

gphn <- 268566 
subg <- make_ego_graph(g, order = 1, nodes = V(g)[names(V(g)) == "Gphn"])[[1]]

# Partition the network.
# The function getProfile 
partition <- getProfile(g,partition_type="ModularityVertexPartition",nsteps=1)

partition_type="ModularityVertexPartition"
nsteps=1


#------------------------------------------------------------------------------
# Finding the shortest path between proteins.
#------------------------------------------------------------------------------

library(getPPIs)

# The mouse interactome.
genes <- unique(c(musInteractome$osEntrezA, musInteractome$osEntrezB))

genes <- compiled_iPSD
netDat <- buildNetwork(musInteractome,genes)
g <- netDat$network

# Find shortest paths between two proteins.
paths <- shortest_paths(g, from="Gphn", to="Nrcam", output = c("vpath", "epath", "both"))


