# Downloading NCBI homology database. And mapping proteins to a species of interst.

# Download and load NCBI homology gene data.
getHomoloGene <- function(hitpredict,taxid,downloads=getwd(),keepdata=FALSE){

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
	# Remove raw data.
	if (keepdata==FALSE) {
		unlink(destfile)
	}
	# Fix column names.
	# Gene ID is a genes organism specific Entrez ID.
	# HID is the genes homology id.
	homologene <- homologene %>% rename_all(list(~ c("HID", "TaxonomyID", "GeneID",
							 "GeneSymbol","ProteinGI","ProteinAccession")))

	# Use HomoloGene to create homology map.
	homology_map <- as.list(homologene$HID)
	names(homology_map) <- homologene$GeneID

	# Map EntrezA to homology ID.
	hidA <- homology_map[hitpredict$EntrezA]
	hidA[seq(1:length(hidA))[unlist(lapply(hidA,is.null))]] <- NA
	hitpredict$HIDA <- unlist(hidA)
	# Map EntrezB to homology ID.
	hidB <- homology_map[hitpredict$EntrezB]
	hidB[seq(1:length(hidB))[unlist(lapply(hidB,is.null))]] <- NA
	hitpredict$HIDB <- unlist(hidB)

	# Subset homologene data.
	homologene <- homologene %>% filter(TaxonomyID==taxid)

	# Taxonomy info.
	source("annotationDBs.R")
	osDB <- unlist(annotationDBs[sapply(annotationDBs,"[",1) == taxid])

	organism <- osDB[3]

	# Get entrez associated with these HIDs.
	colAB <- paste0(organism,c("EntrezA","EntrezB"))

	hitpredict <- dplyr::mutate(hitpredict,colA=homologene$GeneID[match(hitpredict$HIDA,homologene$HID)])
	hitpredict <- dplyr::mutate(hitpredict,colB=homologene$GeneID[match(hitpredict$HIDB,homologene$HID)])
	# Status report.

	is_missing <- is.na(hitpredict$osEntrezA) | is.na(hitpredict$osEntrezB) 
	n <- length(is_missing)
	message(paste(round(100*sum(is_missing)/n,2), "of interactions mapped to homologous genes."))
	# Return data with genes mapped to gene specific entrez.
	return(hitpredict)
}

