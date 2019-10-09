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
	# All genes from your species of interest.
	homologene <- homologene %>% filter(TaxonomyID==taxid)

	# Taxonomy info.
	source("annotationDBs.R")
	osDB <- unlist(annotationDBs[sapply(annotationDBs,"[",1) == taxid])
	names(osDB) <- sapply(strsplit(names(osDB),"\\."),"[",2)
	organism <- osDB["alias"]

	# Get entrez associated with these HIDs.
	hitpredict <- dplyr::mutate(hitpredict,osEntrezA=homologene$GeneID[match(hitpredict$HIDA,homologene$HID)])
	hitpredict <- dplyr::mutate(hitpredict,osEntrezB=homologene$GeneID[match(hitpredict$HIDB,homologene$HID)])

	# Status report.
	is_missing <- is.na(hitpredict$osEntrezA) | is.na(hitpredict$osEntrezB) 
	n <- length(is_missing)
	message(paste0(round(100*sum(is_missing)/n,2), "% of interactions were mapped to homologous genes in ", organism,"."))

	# Remove rows with NA.
	hitpredict <- dplyr::filter(hitpredict, is_missing)
	nppis <- format(dim(hitpredict)[1],nsmal=0,big.mark=",")
	ngenes <- format(length(unique(c(hitpredict$osEntrezA,hitpredict$osEntrezB))),
			 nsmall=0,big.mark=",")
	message(paste0(nppis," protein-protein interactions among ", ngenes,
		      " genes were compiled from ", organism,"."))

	# Return data with genes mapped to gene specific entrez.
	return(hitpredict)
}

