#' buildNetwork
#'
#' Given a HitPredict data.frame (PPIs) and a vector of entrez IDs, build a 
#' protein-protein interaction graph.
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
#' buildNetwork()
#'
buildNetwork <- function(hitpredict, mygenes, taxid = 10090) {
  # Fixme: Only works for mouse!!!
  # Imports.
  suppressPackageStartupMessages({
    require(getPPIs)
    require(dplyr)
    require(org.Mm.eg.db)
    require(igraph)
  })
  # Load proteins of interest.
  if (is.character(mygenes)) {
    # if character, then read file into R.
    myprots <- data.table::fread(mygenes)$Entrez
    idy <- grep("entrez", tolower(colnames(myprots)))
    mygenes <- unlist(myprots[, ..idy])
    names(mygenes) <- NULL
  } else {
    # Use genes passed by user.
    myprots <- mygenes
  }
  # Get interactions among genes of interest.
  ppis <- hitpredict %>% filter(osEntrezA %in% mygenes & osEntrezB %in% mygenes)
  # Create graph.
  sif <- ppis %>% dplyr::select(osEntrezA, osEntrezB)
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
  # Build igraph object.
  g <- graph_from_data_frame(sif, directed = FALSE, vertices = noa)
  g <- simplify(g) # remove any redundant edges.
  nNodes <- length(V(g))
  nEdges <- length(E(g))
  message(paste(nEdges,"identified among","nNodes!"))
  return(g)
}
