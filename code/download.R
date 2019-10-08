#!/usr/bin/env Rscript
# Download HitPredict database.

#------------------------------------------------------------------------------
## Download and compile HitPredict dababase.
#------------------------------------------------------------------------------

download <- function(organism="HitPredict", downloads=getwd()) {
	# Global Imports.
	suppressPackageStartupMessages({
		require(dplyr)
	})
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
		message(paste("Downloading", organism, "PPIs from HitPredict.org."))
		download.file(url,gzfile)
		untar(gzfile, exdir=downloads)
		rawdat <- data.table::fread(myfile,header=TRUE,skip=5)
		unlink(myfile)
	} else {
		message("Using previously downloaded file!")
		untar(gzfile, exdir=downloads)
		rawdat <- data.table::fread(myfile,header=TRUE,skip=5)
		unlink(myfile)
	}
	# Clean up the raw data.
	cleandat <- rawdat %>% rename_all(list(~ gsub(" ","_",.)))
	cleandat$Interactor_A_ID <- gsub("uniprotkb:","",cleandat$Interactor_A_ID)
	cleandat$Interactor_B_ID <- gsub("uniprotkb:","",cleandat$Interactor_B_ID)

	# Organism specific mapping databases:
	annotationDBs <- list(Anopheles_gambia = list(taxid=7165, database="org.Ag.eg.db", alias="anopheles"),
			      Arabidopsis_thaliana = list(taxid=7302, database="org.At.tair.db", alias="arabidopsis"),
			      Bos_taurus = list(taxid=9913, database="org.Bt.eg.db", alias="bovine"),
			      Caenorhabditis_elegans = list(taxid=6239, database="org.Ce.eg.db",alias="worm"),
			      Canis_familiaris = list(taxid=9615,database="org.Cf.eg.db",alias="dog"),
			      Drosophila_melanogaster = list(taxid=7227,database="org.Dm.eg.db",alias="fly"),
			      Danio_rerio = list(taxid=7955, database="org.Dr.eg.db",alias="zebrafish"),
			      Escherichia_coli = list(taxid=83333, database="org.EcK12.eg.db",alias="ecoli-k12"),
			      Escherichia_coliSakai = list(taxid=386585, database="org.EcSakai.eg.db",alias="ecoli-sakai"),
			      Gallus_gallus = list(taxid=9031, database="org.Gg.eg.db",alias="chicken"),
			      Homo_sapiens = list(taxid=9606,database="org.Hs.eg.db",alias="human"),
			      Mus_musculus = list(taxid=10090,database="org.Mm.eg.db",alias="mouse"),
			      Macaca_mulatta = list(taxid=9544,database="org.Mmu.eg.db",alias="macaque"),
			      Plasmodium_falciparum = list(taxid=5833,database="org.Pf.plasmo.db",alias="malaria"),
			      Pan_troglodytes = list(taxid=9598,database="org.Pt.eg.db",alias="chimp"),
			      Rattus_norvegicus = list(taxid=10116,database="org.Rn.eg.db",alias="rat"),
			      Saccharomyces_cerevisiae = list(taxid=4932,database="org.Sc.sgd.db",alias="yeast"),
			      Sus_scrofa = list(taxid=9823,database="org.Ss.eg.db",alias="pig"),
			      Xenopus_laevis= list(taxid=8355,database="org.Xl.eg.db",alias="frog"))

	# Subset HitPredict data, keep interactions from species with mapping databases..
	dbs <- unlist(sapply(annotationDBs,"[",1))	

	data <- data.table::as.data.table({
		cleandat %>% filter(Interactor_A_Taxonomy %in% dbs & Interactor_B_Taxonomy %in% dbs)})

	# Taxonomy ids for remaining organisms.	
	taxids <- unique(c(data$Interactor_A_Taxonomy, data$Interactor_B_Taxonomy))
	names(taxids) <- unlist(sapply(annotationDBs,"[",3))[match(taxids, unlist(sapply(annotationDBs,"[",1)))]

	# Loop to map UniprotIDs to Entrez for all species.
	for (species in taxids){
		# Get uniprot IDs.
		subdat <- data %>% dplyr::filter(Interactor_A_Taxonomy == species & Interactor_B_Taxonomy == species)
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
		entrez <- AnnotationDbi::mapIds(myDB,
					 keys = uniprot, 
					 column = "ENTREZID", 
					 keytype = "UNIPROT", 
					 multiVals="first")
		# Status report.
		percent_mapped <- round(100*(1-sum(is.na(entrez))/length(entrez)),2)
		message(paste(percent_mapped, "of",names(taxids)[taxids==species], 
			"Uniprot IDs were successfully mapped to Entrez IDs."))
		# Create protein identifier map.
		protMap <- as.list(entrez)
		names(protMap) <- uniprot
		# Add to HitPredict data.
		subdat$EntrezA <- unlist(protMap[subdat$Interactor_A_ID])
		subdat$EntrezB <- unlist(protMap[subdat$Interactor_B_ID])
		# Add back to data.
		data <- left_join(data,subdat, by = c("Interactor_A_ID", "Interactor_B_ID", 
						      "Interactor_A_Name", "Interactor_B_Name", 
						      "Interactor_A_Alias", "Interactor_B_Alias", 
						      "Interaction_detection_methods", 
						      "First_author", "Publications", 
						      "Interactor_A_Taxonomy","Interactor_B_Taxonomy", 
						      "Interaction_type", "Source_database", 
						      "Interaction_identifier", "Confidence_score"))
	} # ENDS loop
	return(data)
}


foo <- download()
	



