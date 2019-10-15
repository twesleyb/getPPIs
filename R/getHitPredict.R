#------------------------------------------------------------------------------
#' getHitPredict
#'
#' Function that facilitates downloading the HitPredict database.
#'
#' @param 
#'
#' @return 
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references 
#'
#' @keywords 
#'
#' @export 
#'
#' @examples 

getHitPredict <- function(organism="HitPredict", downloads=getwd(), keepdata=FALSE) {
	# Imports.
	suppressPackageStartupMessages({
		require(dplyr)})
	# Parse HitPredict Downloads page.
	url <- "http://www.hitpredict.org/download.html"
	pg <- xml2::read_html(url)
	nodes <- rvest::html_nodes(pg,"a")
	hrefs <- rvest::html_attr(nodes, "href")
	hrefs <- sub("./","", hrefs[grepl("MITAB-2.5.tgz",hrefs)])
	namen <- sapply(strsplit(basename(hrefs),"_interactions"),"[",1)
	hrefs <- as.list(hrefs)
	names(hrefs) <- namen
	# Download HitPredict database as MITAB-2.5 format.
	url <- file.path("http://www.hitpredict.org", hrefs[[organism]])
	gzfile <- file.path(downloads, basename(url))
	myfile <- tools::file_path_sans_ext(gzfile)
	if (!file.exists(gzfile)) {
		message(paste("Downloading", organism, "PPIs from HitPredict.org..."))
		download.file(url,gzfile)
		untar(gzfile, exdir=downloads)
		rawdat <- data.table::fread(myfile,header=TRUE,skip=5)
	} else {
		message("Using previously downloaded HitPredict data!")
		untar(gzfile, exdir=downloads)
		rawdat <- data.table::fread(myfile,header=TRUE,skip=5)
	}
	if (keepdata==FALSE) {
		unlink(gzfile)
		unlink(myfile)
	}
	# Clean up the raw data.
	cleandat <- rawdat %>% rename_all(list(~ gsub(" ","_",.)))
	cleandat$Interactor_A_ID <- gsub("uniprotkb:","",cleandat$Interactor_A_ID)
	cleandat$Interactor_B_ID <- gsub("uniprotkb:","",cleandat$Interactor_B_ID)
	# Define organism specific mapping databases:
	source("annotationDBs.R")
	# Subset HitPredict data, keep interactions from species with mapping databases...
	dbs <- unlist(sapply(annotationDBs,"[",1))	
	data <- data.table::as.data.table({
		cleandat %>% filter(Interactor_A_Taxonomy %in% dbs & Interactor_B_Taxonomy %in% dbs)})
	data$EntrezA <- "NA"
	data$EntrezB <- "NA"
        # The taxonnomy ids for remaining organisms:
	taxids <- unique(c(data$Interactor_A_Taxonomy, data$Interactor_B_Taxonomy))
	names(taxids) <- unlist(sapply(annotationDBs,"[",3))[match(taxids, unlist(sapply(annotationDBs,"[",1)))]
	# Loop to map UniprotIDs to Entrez for all species (taxids).
	message("Mapping Uniprot IDs to Entrez IDs...")
	for (species in taxids){
		# Get uniprot IDs.
		subdat <- data %>% dplyr::filter(Interactor_A_Taxonomy == species & Interactor_B_Taxonomy == species)
		filtdat <- data %>% dplyr::filter(Interactor_A_Taxonomy != species & Interactor_B_Taxonomy != species)
		uniprot <- subdat %>% dplyr::select(Interactor_A_ID, Interactor_B_ID) %>% 
			stack() %>% pull(values)
		# Get organism specific mapping database.
		orgDB <- unlist(annotationDBs[sapply(annotationDBs,"[",1)==species])
		names(orgDB) <- sapply(strsplit(names(orgDB),"\\."),"[",2)
		suppressPackageStartupMessages({
			eval(parse(text=paste0("require(",orgDB[["database"]],",quietly=TRUE)")))
		})
		myDB <- eval(parse(text=orgDB[["database"]]))
		# Perform Uniprot 3 Entrez mapping
		suppressMessages({
			entrez <- AnnotationDbi::mapIds(myDB,
					 keys = uniprot, 
					 column = "ENTREZID", 
					 keytype = "UNIPROT", 
					 multiVals="first")
		})
		# Status report.
		percent_mapped <- round(100*(1-sum(is.na(entrez))/length(entrez)),2)
		message(paste("   ...", percent_mapped, "% of", names(taxids)[taxids==species], 
			"Uniprot IDs were successfully mapped to Entrez IDs."))
		# Create protein identifier map.
		protMap <- as.list(entrez)
		names(protMap) <- uniprot
		# Add to HitPredict data.
		subdat$EntrezA <- unlist(protMap[subdat$Interactor_A_ID])
		subdat$EntrezB <- unlist(protMap[subdat$Interactor_B_ID])
		data <- full_join(filtdat,subdat, by = c("Interactor_A_ID", "Interactor_B_ID", 
						      "Interactor_A_Name", "Interactor_B_Name", 
						      "Interactor_A_Alias", "Interactor_B_Alias", 
						      "Interaction_detection_methods", 
						      "First_author", "Publications", 
						      "Interactor_A_Taxonomy","Interactor_B_Taxonomy", 
						      "Interaction_type", "Source_database", 
						      "Interaction_identifier", "Confidence_score",
						      "EntrezA","EntrezB"))
	}
	# Progress Report.
	N <- dim(data)[1]
	n_missing <- sum(is.na(data$EntrezA) | is.na(data$EntrezB))
	percent_missing <- round(100*n_missing/N,2)
	n_missing <- format(as.numeric(n_missing), nsmall=0, big.mark=",")
	message(paste0(n_missing, " (", percent_missing, "%) ", 
		      "of all interactions (rows) are incompletely mapped to Entrez IDs and will be discarded."))
	# Remove unmapped ids.
	data <- dplyr::filter(data, !is.na(EntrezA) & !is.na(EntrezB))
	n <- format(as.numeric(dim(data)[1]), nsmall=0, big.mark=",")
	message(paste(n,"Protein-protein interactions from",
		     length(taxids), "species were compiled from the HitPredict",
		     "database and successfully mapped to Entrez IDs."))
	# Return HitPredict data with Entrez IDs.
	return(data)
}
#------------------------------------------------------------------------------
