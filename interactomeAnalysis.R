# Analysis of the mouse interactome.

# Load the data.
library(getPPIs)
data(musInteractome)

# Build a network.
mygenes <- unique(c(musInteractome$osEntrezA,musInteractome$osEntrezB))
g <- buildNetwork(musInteractome,mygenes)


#Error in graph <- from <- data <- frame(sif, directed = FALSE, vertices = noa) :
#some vertex names in edge list are not listed in vertex data frame


# Node names.
nodes <- names(V(g))
names(nodes) <- vertex_attr(g,"symbol")

# Partition the network.
partition <- getProfile(g,partition_type="ModularityVertexPartition",nsteps=1)

#------------------------------------------------------------------------------

