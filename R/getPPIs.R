#' getPPIs
#'
#' Wrapper around getHitPredict, getHomoloGene, and getMethods.
#' Provided a list of genes, get PPIs amongst these proteins.
#'
#' @param organism (character) organism to be downloaded from HitPredict.
#' One of c("HitPredict",...)
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
getPPIs <- function(organism = "HitPredict",
                    mygenes,
                    taxid,
                    methods = "all",
                    cutoff = 0,
                    downloads = getwd(),
                    keepdata = FALSE,
                    saveppis = FALSE) {
  # Imports.
  suppressPackageStartupMessage({
    require(getPPIs)
    require(dplyr)
    require(org.Mm.eg.db)
  })
  # Download HitPredict database.
  hitpredict <- getHitPredict(organism, downloads, keepdata)
  # Map genes to homologous mouse genes.
  hitpredict <- getHomoloGene(hitpredict, taxid, downloads, keepdata)
  # Annotate hitpredict data with method names.
  hitpredict <- getMethods(hitpredict, methods, cutoff)
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
  # Write to file.
  if (saveppis == TRUE) {
    message("Saving compiled protein-protein interactions to file!")
    data.table::fwrite(ppis, file.path(paste0(file, "PPIs.csv")))
  }
  # Return all ppis among proteins of interest.
  return(ppis)
}
