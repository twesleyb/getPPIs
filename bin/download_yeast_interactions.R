#!/usr/bin/env Rscript

# INPUT in root/downloads: 
# NOTE: Data were downloaded as tsv from YeastMine on 08/03/2020. 
#       See xml queries in root/docs.
# NOTE: Column headers where edited by hand in excel and saved as txt.
data_files <- c("Yeast_Hs_Homologs_08032020.txt", # Yeast gene human homologs.
		"Yeast_Interactions_08032020.txt") # Yeast interactions.

# Load renv.
here <- getwd()
root <- dirname(here)
renv::load(root)

# Load the data.
downdir <- file.path(root,"downloads")
homologs <- data.table::fread(file.path(downdir,data_files[1]))
interactions <- data.table::fread(file.path(downdir,data_files[2]))

# Save the data as rda.
save(homologs,file="yeast_homologs.rda",version=2)
save(interactions,file="yeast_interactions.rda",version=2)
