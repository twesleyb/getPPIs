# Analysis of the mouse interactome.


# Load the data.
library(getPPIs)
data(musInteractome)

# Build a network.
mygenes <- unique(c(musInteractome$osEntrezA,musInteractome$osEntrezB))
g <- buildNetwork(musInteractome,mygenes)

# Node names.
nodes <- names(V(g))
names(nodes) <- vertex_attr(g,"symbol")

# Partition the network.
getProfile(g,partition_type="ModularityVertexPartition",nsteps=1)

#------------------------------------------------------------------------------

