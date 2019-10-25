#' mapIDs
#'
#' Mapping gene identifiers.
#'
#' @param identifiers
#'
#' @param from - One of c("ACCNUM", "ALIAS", "ENSEMBL", "ENSEMBLPROT",
#' "ENSEMBLTRANS", "ENTREZID", "ENZYME", "EVIDENCE", "EVIDENCEALL", "GENENAME",
#' "GO", "GOALL", "IPI", "MGI", "ONTOLOGY", "ONTOLOGYALL", "PATH",
#' "PFAM", "PMID", "PROSITE", "REFSEQ", "SYMBOL", "UNIGENE", "UNIPROT")
#'
#' @param to - see from
#'
#' @param species
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
#' @examples
#' mapIDs(mygenes, from = "symbols", to = "entrez", species = "mouse")
mapIDs <- function(identifiers, from, to, species, ...) {
  # Wrapper around AnnotationDbi::mapIds()
  # Imports.
  require(getPPIs)
  # Get organism specific mapping database.
  annotationDBs <- mappingDBs()
  orgDB <- unlist(annotationDBs[sapply(annotationDBs, "[", 3) == tolower(species)])
  names(orgDB) <- sapply(strsplit(names(orgDB), "\\."), "[", 2)
  suppressPackageStartupMessages({
    eval(parse(text = paste0("require(", orgDB[["database"]], ",quietly=TRUE)")))
  })
  osDB <- eval(parse(text = orgDB[["database"]]))
  # Map gene identifiers. Suppress messages.
  suppressMessages({
    output <- AnnotationDbi::mapIds(osDB,
      keys = as.character(identifiers),
      column = columns(osDB)[grep(toupper(to), columns(osDB))],
      keytype = columns(osDB)[grep(toupper(from), columns(osDB))],
      multiVals = "first",
    )
  })
  # Check that all nodes (entrez) are mapped to gene symbols.
  not_mapped <- is.na(output)

  if (sum(is.na(output)) != 0) {
    message(paste0("Warning: Unable to map ", sum(not_mapped)," ",species," ",
		   from,"(s)"," to ",to," identifiers!"))
  }
  names(output) <- identifiers
  return(output)
}
