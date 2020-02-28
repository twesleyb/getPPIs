#' getPPIs
#'
#' Wrapper around getHitPredict, getHomoloGene, and getInteractionMethods.
#'
#' @param dataset (character) dataset to be downloaded from HitPredict.
#' One of c("HitPredict",...)
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
getPPIs <- function(dataset, species) {
  # Download HitPredict database.
  hitpredict <- getHitPredict(dataset)
  # Map genes to homologous mouse genes.
  taxid <- getTaxid(species)
  hitpredict <- getHomoloGene(hitpredict, taxid)
  # Annotate hitpredict data with method names.
  hitpredict <- getInteractionMethods(hitpredict)
  # Status.
  message("Complete!")
  return(hitpredict)
}
