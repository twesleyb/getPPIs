#!/usr/bin/env Rscript
#-------------------------------------------------------------------------------

# Load functions.
source("getHitPredict.R")
source("getHomoloGene.R")

# Directories.
here <- getwd()
rootdir <- dirname(here)
downloads <- file.path(rootdir,"downloads")

# Download HitPredict database.
hitpredict <- getHitPredict("HitPredict", downloads)

# Download Homology data.
getHomoloGene()

# Create homology map for human, mouse, and rat.
homology_map <- as.list(homologene$HID)
names(homology_map) <- homologene$GeneID

# Map EntrezB to homology ID.
hidA <- homology_map[data$EntrezA]
hidA[seq(1:length(hidA))[unlist(lapply(hidA,is.null))]] <- NA
data$HIDA <- unlist(hidA)
# Map EntrezA to homology ID.
hidB <- homoMap[data$EntrezB]
hidB[seq(1:length(hidB))[unlist(lapply(hidB,is.null))]] <- NA
data$HIDB <- unlist(hidB)

#------------------------------------------------------------------------------
# Map HIDs to mouse homologous gene.
#------------------------------------------------------------------------------

# Download mouse homology data from MGI.
url <- "http://www.informatics.jax.org/downloads/reports/HGNC_homologene.rpt"
myfile <- file.path(downloads, "HGNC_homologene.rpt")
if (!file.exists(myfile)) {
	message("Downloading Mouse Homology data from MGI...")
	download.file(url, myfile)
	data <- data.table::fread(myfile, header=FALSE)
} else {
	message("file already exists!")
	# Use read.delim, bc fread seems to not get it right.
	data <- data.table::fread(myfile, header=FALSE)
}

# Clean up the data...
data <- data %>% dplyr::select(-ncol(data)) %>%
	rename_all(list(~ names(read.delim(myfile,nrow=1))))



