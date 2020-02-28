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
getPPIs <- function(dataset="all", species) {
  # Download HitPredict database.
  hitpredict <- getHitPredict(dataset)
  # Map genes to homologous mouse genes.
  hitpredict$osEntrezA <- getHomologs(hitpredict$EntrezA,species)
  hitpredict$osEntrezB <- getHomologs(hitpredict$EntrezB,species)
  # Remove missing ids.
  is_missing <- is.na(hitpredict$osEntrezA) | is.na(hitpredict$osEntrezB)
  hitpredict <- subset(hitpredict,!is_missing)
  # Status report.
  nPPIs <- formatC(nrow(hitpredict),big.mark=",")
  message(paste(nPPIs,"Protein-protein interactions were successfully",
	       	"mapped to", species,"Entrez IDs."))
  # Annotate hitpredict data with method names.
  hitpredict <- getInteractionMethods(hitpredict)
  # Status.
  message("Complete!")
  return(hitpredict)
}
