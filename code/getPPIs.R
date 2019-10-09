#!/usr/bin/env Rscript

# Load functions.
source("getHitPredict.R")
source("getHomoloGene.R")

# Directories.
here <- getwd()
rootdir <- dirname(here)
downloads <- file.path(rootdir,"downloads")

# Download HitPredict database.
hitpredict <- getHitPredict("HitPredict")

# Map genes to homologous gene of a species of interest.
ppis <- getHomoloGene(hitpredict, taxid = 10090) # Mouse

# Create a network.
net <- getNetwork(mygenes,ppis)
