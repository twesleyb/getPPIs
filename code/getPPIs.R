#!/usr/bin/env Rscript

# Wrapper around getHitPredict, getHomoloGene, and getMethods.
# Provided a list of genes, get PPIs.

# Directories.
here <- getwd()
rootdir <- dirname(here)
downloads <- file.path(rootdir,"downloads")
datadir <- file.path(rootdir,"data")
funcdir <- file.path(rootdir,"R")

# Load functions.
#devtools::install_github("twesleyb/getPPIs") # installation
library(getPPIs)

# Download HitPredict database.
hitpredict <- getHitPredict("HitPredict", downloads, keepdata=TRUE)

# Map genes to homologous mouse genes.
hitpredict <- getHomoloGene(hitpredict, taxid = 10090, downloads, keepdata=TRUE)

# Annotate hitpredict data with method names.
# Fix to automate download of mi.owl or save as Rdata and import.
path2mi <- file.path(downloads,"mi.owl")
hitpredict <- getMethods(hitpredict,path2mi, methods = "all", cutoff=0)

# Write to file.
#data.table::fwrite(hitpredict,file.path(datadir,"musInteractome.csv"))

# Load compiled iPSD proteome.
# Create a network.
myfile <- file.path(prots,"compiled_iPSD.txt")
myfile <- file.path(prots,"ePSD.txt")
mygenes <- data.table::fread(myfile)$Entrez

# Get interactions among genes of interest.
ppis <- hitpredict %>% filter(osEntrezA %in% mygenes & osEntrezB %in% mygenes)

# Write to file.
data.table::fwrite(ppis,"ePSD_PPIs.csv")


# Simple interaction file:
sif <- ppis %>% dplyr::select(osEntrezA, osEntrezB)

# Create node attribute data.table.
nodes <- unique(c(sif$osEntrezA,sif$osEntrezB))
symbols <- AnnotationDbi::mapIds(org.Mm.eg.db,
				 keys=as.character(nodes),
				 column="SYMBOL",
				 keytype="ENTREZID",
				 multiVals="first") 
noa <- data.table::data.table(node=nodes,symbol=symbols)

# Check.
if (sum(is.na(symbols)) == 0) { 
	message("All node Entrez IDs mapped to gene symbols!") 
}


# Save data.
#data.table::fwrite(noa,file.path(datadir,"noa.csv"))
#data.table::fwrite(sif,file.path(datadir,"sif.csv"))

# Create igraph.
g <- igraph::graph_from_data_frame(sif, directed = FALSE, vertices = noa)

# Simplify to remove redundant edges.
g <- simplify(g,remove.loops=TRUE)

# Basic graph properties.
length(V(g))
length(E(g))

# Evaluate topology...
# Need to fix this. Examination of Cytoscape stats suggest that ther is scale free topology.
library(TBmiscr)
k <- degree(g, loops=FALSE)
plot <- ggplotScaleFreePlot(k)
plot

#----------------------------------------------------
# Community detection!
library(leiden)

# Single partition.
# Why do some methods not work?
partition <- leiden(g, partition_type = "RBConfigurationVertexPartition", seed = 1, n_iterations = -1)

# Resolution profile.
profile <- list()
res <- seq(0,1,0.01)
for (i in 1:length(res)){
	print(i)
	r <- res[i]
	profile[[i]] <- leiden(g,partition_type = "RBConfigurationVertexPartition", resolution_parameter = r, seed = 1, n_iterations=-1)
}

# Using lapply.
profile <- lapply(as.list(seq(0,1,0.01)), function(x) leiden(g, partition_type = "RBConfigurationVertexPartition",
		  resolution_parameter = x, seed = 1, n_iterations = -1))

# Number of clusters at each resolution.
k = unlist(lapply(profile,function(x) length(unique(x))))
names(k) <- paste0("r=",seq(0,1,0.01))

# Modularity of each partition.
Q <- unlist(lapply(profile, function(x) modularity(g, membership=x)))
names(Q) <- paste0("r=",seq(0,1,0.01))

# Best params.
best_resolution <- names(Q[Q==max(Q)])
Qmax <- Q[Q==max(Q)]
koptimum <- k[best_resolution]

best_resolution
Qmax
koptimum
