#' getPPIs
#'
#' Wrapper around getHitPredict, getHomoloGene, and getMethods.
#'
#' @param organism (character) organism to be downloaded from HitPredict.
#' One of c("HitPredict",...)
#'
#' @param taxid (integer) taxonomic identifier for organism of interest.
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
getPPIs <- function(organism, taxid) {
  # Download HitPredict database.
  hitpredict <- getHitPredict(organism)
  # Map genes to homologous mouse genes.
  hitpredict <- getHomoloGene(hitpredict, taxid)
  # Annotate hitpredict data with method names.
  hitpredict <- getMethods(hitpredict)
  # Status.
  message("Complete!")
  return(hitpredict)
}
