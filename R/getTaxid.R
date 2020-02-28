#' getTaxid
#
#' taxonomic id given a species alias
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @export
#'
getTaxid <- function(species){
	annotationDBs <- mappingDBs()
	orgDB <- unlist(annotationDBs[sapply(annotationDBs, 
					     "[", 3) == tolower(species)])
	species <- orgDB$alias
	return(species)
}
