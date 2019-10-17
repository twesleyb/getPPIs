#' getPPIs
#'
#' Wrapper around getHitPredict, getHomoloGene, and getMethods.
#' Provided a list of genes, get PPIs amongst these proteins.
#'
#' @param none
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
getPPIs <- function(organism = "HitPredict", mygenes, taxid = 10090, prots, methods = "all", cutoff = 0,
                    file = "", downloads = getwd(), keepdata = FALSE) {
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
  myprots <- data.table::fread(mygenes)$Entrez
  # Get interactions among genes of interest.
  ppis <- hitpredict %>% filter(osEntrezA %in% mygenes & osEntrezB %in% mygenes)
  # Write to file.
  data.table::fwrite(ppis, file.path(paste0(file, "PPIs.csv")))
  # Simple interaction file:
  sif <- ppis %>% dplyr::select(osEntrezA, osEntrezB)
  # Create node attribute data.table.
  nodes <- unique(c(sif$osEntrezA, sif$osEntrezB))
  symbols <- AnnotationDbi::mapIds(org.Mm.eg.db,
    keys = as.character(nodes),
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  noa <- data.table::data.table(node = nodes, symbol = symbols)
  # Check.
  if (sum(is.na(symbols)) == 0) {
    message("All node Entrez IDs mapped to gene symbols!")
  }
  # Save data.
  data.table::fwrite(noa, file = file.path(rootdir, "temp", "noa.csv"))
  data.table::fwrite(sif, file = file.path(rootdir, "temp", "sif.csv"))
  # Return all ppis among proteins of interest.
  return(ppis)
}
