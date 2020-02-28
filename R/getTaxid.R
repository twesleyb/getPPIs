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
	orgDB <- annotationDBs[sapply(annotationDBs,"[[", 3) == tolower(species)]
	idx <- names(which(tolower(species)==sapply(annotationDBs,"[[",3)))
	if (length(idx)==0) { stop("Please provide a valid species alias.") }
	taxid <- annotationDBs[[idx]][["taxid"]]
	return(taxid)
}
