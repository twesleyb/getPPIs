#' getTaxid
getTaxid <- function(species){
	annotationDBs <- mappingDBs()
	orgDB <- unlist(annotationDBs[sapply(annotationDBs, 
					     "[", 3) == tolower(species)])
	species <- orgDB$alias
	return(species)
}
