#!/usr/bin/env Rscript

## Build a human interactome using the getPPIs package.

# Load the library.
library(getPPIs)

# Build the human interactome.
musInteractome <- getPPIs("HitPredict", species="mouse")

# Summarize total number of interactions.
genes <- unique(c(musInteractome$EntrezA,musInteractome$EntrezB))
nGenes <- formatC(length(genes),big.mark=",")
nPPIs <- formatC(nrow(musInteractome),big.mark=",")
message(paste("Total number of Genes:",nGenes))
message(paste("Total number of PPIS:",nPPIs))

# Number of publications.
publications <- unique(unlist(strsplit(musInteractome$Publications,"\\|")))
nPubs <- formatC(length(publications),big.mark=",")

# Number of species.
nSpecies <- length(unique(musInteractome$Interactor_A_Taxonomy))

# Number of interactions from each species.
df <- t(table(musInteractome$Interactor_A_Taxonomy))

# Map taxid to species alias.
dbs <- mappingDBs()
taxids <- sapply(dbs,"[",3)
names(taxids) <- sapply(dbs,"[",1)

# Summary.
colnames(df) <- taxids[colnames(df)]
df <- apply(df,1,function(x) formatC(x,big.mark=","))
message("Number of interactions from each species:")
knitr::kable(t(df))

# Table summarizing number of interactions from each database.
df <- t(table(musInteractome$Source_database))
df <- apply(df,1,function(x) formatC(x,big.mark=","))
knitr::kable(t(df))

# Save as rda object.
save(musInteractome,file="musInteractome.rda")

# In my experience, when working with mouse proteins, it is useful to
# only consider interactions identified in mouse, human, and/or rat.
orgs <- c(10090,9606,10116)
idx <- musInteractome$Interactor_A_Taxonomy %in% orgs
ppis <- subset(musInteractome, idx)
