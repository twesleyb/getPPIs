#!/usr/bin/env Rscript

# Load functions.
source("getHitPredict.R")
source("getHomoloGene.R")

# Directories.
here <- getwd()
rootdir <- dirname(here)
downloads <- file.path(rootdir,"downloads")
datadir <- file.path(rootdir,"data")

# Download HitPredict database.
hitpredict <- getHitPredict("HitPredict", downloads, keepdata=TRUE)

# Map genes to homologous mouse genes.
hitpredict <- getHomoloGene(hitpredict, taxid = 10090, downloads, keepdata=TRUE)

# Write to file.
data.table::fwrite(hitpredict,file.path(datadir,"musInteractome.csv"))

quit()

# Load compiled iPSD proteome.
# Create a network.
mygenes <-
net <- getNetwork(mygenes,ppis)
