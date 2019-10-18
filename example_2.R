#------------------------------------------------------------------------------
## Using the pre-built mouse Interactome.
#------------------------------------------------------------------------------

devtools::install_github("twesleyb/getPPIs")

library(getPPIs)

# Load the mouse interactome.
data(musInteractome)
data(compiled_iPSD)

# Build Synaptosome PPI graph.
g <- buildNetwork(musInteractome, compiled_iPSD)

# Community detection with the Leiden algorithm.
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



