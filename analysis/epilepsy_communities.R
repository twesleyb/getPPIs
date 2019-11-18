# Build networks of proteins associated with epilepsy genes.

library(getPPIs)

# Load epilepsy-associated genes and mouse interactome.
data(epilepsy_genes)
data(musInteractome)

# Get mouse homologs to epilepsy genes.
genes <- getHomologs(epilepsy_genes$all_epilepsy_genes,10090)
seeds <- genes[!is.na(genes)]

# Build a mouse interactome.
all_genes <- unique(c(musInteractome$osEntrezA,musInteractome$osEntrezB))
graph <- buildNetwork(musInteractome,all_genes,10090)

# Add community nodes, proteins that form >2 connections to epilepsy genes.
subg <- getCommunity(graph,seeds)

# Community detection.


