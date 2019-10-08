# Downloading NCBI homology database.

# Download and load NCBI homology gene data.
getHomoloGene <- function(){
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
							 "GeneSymbol","ProteinGI","ProteinAccession")))
	return(homologene)
}

