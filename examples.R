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
data(compiled_iPSD)

# Build some graphs.
ipsd <- buildNetwork(hitpredict=musInteractome, mygenes=compiled_iPSD, taxid=10090)
wrp <- buildNetwork(hitpredict=musInteractome, mygenes=Wrp, taxid=10090)
epsd <- buildNetwork(hitpredict=musInteractome, mygenes=ePSD, taxid=10090)


plot(ipsd$network)

plot(wrp$network)

plot(epsd$network) 

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


