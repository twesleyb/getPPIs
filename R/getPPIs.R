#' getPPIs
#'
#' Wrapper around getHitPredict, getHomoloGene, and getMethods.
#' Provided a list of genes, get PPIs amongst these proteins.
#'
#' @param organism (character) organism to be downloaded from HitPredict.
#' One of c("HitPredict",...)
#'
#' @param taxid (integer) taxonomic identifier for organism of interest.
#'
#' @param mygenes (character or vector of integers) a character specifying the
#' file containing Entrez IDs for your genes of interest. Should contain a column,
#' Entrez with your genes. Alternatively, can be a vector of integers, Entrez IDs.
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
getPPIs <- function(organism, taxid, mygenes) {
  # Imports.
  suppressPackageStartupMessages({
    require(getPPIs)
    require(dplyr)
  })
  # Download HitPredict database.
  hitpredict <- getHitPredict(organism)
  # Map genes to homologous mouse genes.
  hitpredict <- getHomoloGene(hitpredict, taxid)
  # Annotate hitpredict data with method names.
  hitpredict <- getMethods(hitpredict, methods = "all", cutoff)
  # Load proteins of interest.
  if (is.character(mygenes)) {
    # if character, then read file into R.
    myprots <- data.table::fread(mygenes)$Entrez
  } else {
    # Use genes passed by user.
    myprots <- mygenes
  }
  # Get interactions among genes of interest.
  ppis <- hitpredict %>% filter(osEntrezA %in% mygenes & osEntrezB %in% mygenes)
  # Return all ppis among proteins of interest.
  return(ppis)
}
