#!/usr/bin/env Rscript

# Load functions.
source("download.R")

# Directories.
here <- getwd()
rootdir <- dirname(here)
downloads <- file.path(rootdir,"downloads")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Download HitPredict database.
hitpredict <- download("HitPredict", downloads)

# Download and load NCBI homology gene data.
url <- "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"
destfile <- file.path(downloads, "homologene.data")
if (!file.exists(destfile)) {
	message("Downloading NCBI HomoloGene data...")
	download.file(url, destfile)
	homologene <- data.table::fread(destfile, header = FALSE)
} else {
	message("Loading NCBI HomoloGene data from file!")
	homologene <- data.table::fread(destfile, header = FALSE)
}

# Fix column names.
# Gene ID is a genes organism specific Entrez ID.
# HID is the genes homology id.
homologene <- homologene %>% rename_all(list(~ c("HID", "TaxonomyID", "GeneID",
						 "GeneSymbol", "ProteinGI", "ProteinAccession")))

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



