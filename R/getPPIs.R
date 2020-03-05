#' getPPIs
#'
#' Wrapper around getHitPredict, getHomologs, and getInteractionMethods.
#'
#' @param dataset (character) dataset to be downloaded from HitPredict.
#'
#' @param species (character) species alias for organism of interest.
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' getPPIs()
getPPIs <- function(dataset="all", species, quiet = FALSE) {
  # Download HitPredict database.
  hitpredict <- getHitPredict(dataset,quiet)
  # Map genes to homologous mouse genes.
  hitpredict$osEntrezA <- getHomologs(hitpredict$EntrezA,species,quiet)
  hitpredict$osEntrezB <- getHomologs(hitpredict$EntrezB,species,quiet)
  # Remove missing ids.
  is_missing <- is.na(hitpredict$osEntrezA) | is.na(hitpredict$osEntrezB)
  hitpredict <- subset(hitpredict,!is_missing)
  # Status report.
  nPPIs <- formatC(nrow(hitpredict),big.mark=",")
  if (!quiet) {
	  message(paste(nPPIs,"Protein-protein interactions were successfully",
			"mapped to", species,"Entrez IDs."))
  }
  # Annotate hitpredict data with method names.
  hitpredict <- getInteractionMethods(hitpredict)
  # Status.
  if (!quiet) { message("Complete!") }
  return(hitpredict)
}
