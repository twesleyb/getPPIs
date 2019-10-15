#' getPPIs(taxid, goi, file="", downloads, path2mi)
#'
#' Wrapper around getHitPredict, getHomoloGene, and getMethods.
#' Provided a list of genes, get PPIs.
#'
#' @param none
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @export
#'
#' @examples
#' getPPIs()

# Directories.
here <- getwd()
rootdir <- here
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

# Save
#musInteractome <- hitpredict
#save(musInteractome,file=file.path(datadir,"musInteractome.RData"))

# Load compiled mouse interactome.
data(musInteractome)
hitpredict <- musInteractome

# Load proteins of interest.
myfile <- file.path(rootdir,"temp","paps.txt")
mygenes <- data.table::fread(myfile)$Entrez

# Get interactions among genes of interest.
library(dplyr)
ppis <- hitpredict %>% filter(osEntrezA %in% mygenes & osEntrezB %in% mygenes)

# Write to file.
data.table::fwrite(ppis,file.path(paste0(file,"PPIs.csv")))

# Simple interaction file:
sif <- ppis %>% dplyr::select(osEntrezA, osEntrezB)

# Create node attribute data.table.
nodes <- unique(c(sif$osEntrezA,sif$osEntrezB))
library(org.Mm.eg.db)
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
data.table::fwrite(noa,file=file.path(rootdir,"temp","noa.csv"))
data.table::fwrite(sif,file=file.path(rootdir,"temp","sif.csv"))

# Create igraph.
library(igraph)
g <- graph_from_data_frame(sif, directed = FALSE, vertices = noa)
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

# Save partition.
library(data.table)

fwrite(data.table(Node=names(V(g)),
		  Cluster = partition), 
