#' getIDs
#'
#' Mapping gene identifiers.
#' An easy-to-use and robust wrapper around AnnotationDbi::mapIds().
#'
#' @param identifiers - input gene identifiers
#'
#' @param from - input identifier type, one of (case insensitive):
#' `c("ACCNUM", "ALIAS", "ENSEMBL", "ENSEMBLPROT",
#' "ENSEMBLTRANS", "ENTREZID", "ENZYME", "EVIDENCE", "EVIDENCEALL", "GENENAME",
#' "GO", "GOALL", "IPI", "MGI", "ONTOLOGY", "ONTOLOGYALL", "PATH",
#' "PFAM", "PMID", "PROSITE", "REFSEQ", "SYMBOL", "UNIGENE", "UNIPROT")`
#'
#' @param to - output identifier type, see `from`.
#'
#' @param species - organism identifier for input genes.
#'
#' @return output gene identifiers
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @importFrom AnnotationDbi mapIds
#'
#' @importFrom data.table data.table
#'
#' @export
#'
#' @examples
#' getIDs(mygenes, from = "symbol", to = "entrez", species = "mouse")
getIDs <- function(identifiers, from, to, species = NULL, taxid = NULL,
                   quiet = TRUE, multiVals = "first", ...) {
  # Wrapper around AnnotationDbi::mapIds()
  # Check input identifiers.
  if (sum(is.na(identifiers))) {
    message("Warning: missing values (NA) detected in input identifiers.")
  }
  # Get organism specific mapping database.
  annotationDBs <- mappingDBs()
  if (!is.null(taxid)) {
    orgDB <- unlist(annotationDBs[sapply(annotationDBs, "[", 1) == taxid])
  } else if (!is.null(species)) {
    orgDB <- unlist(annotationDBs[sapply(annotationDBs, "[", 3) == tolower(species)])
  } else {
    stop("Please provide a species or taxid for gene identifiers.")
  }
  names(orgDB) <- sapply(strsplit(names(orgDB), "\\."), "[", 2)
  suppressPackageStartupMessages({
    eval(parse(text = paste0("require(", orgDB[["database"]], ",quietly=TRUE)")))
  })
  osDB <- eval(parse(text = orgDB[["database"]]))
  # Get input type (from) and output type (to).
  colIDto <- grep(toupper(to), columns(osDB))
  colIDfrom <- grep(toupper(from), columns(osDB))
  # Check that from and to map to a single column.
  keys <- keytypes(osDB)
  msg <- paste0(
    "Please provide one of the following",
    "(case insensitive): ", paste(keys, collapse = ", "), "."
  )
  if (length(colIDto) > 1) {
    stop(paste("Multiple matching keys (to).\n", msg))
  }
  if (length(colIDfrom) > 1) {
    stop(paste("Multiple matching keys (from).\n", msg))
  }
  # Check MGI format if input is MGI.
  if (columns(osDB)[colIDfrom] == "MGI") {
    if (!any(grepl("MGI:", identifiers))) {
      stop("Please provide MGI identifiers as MGI:ID")
    }
    identifiers <- paste0(
      "MGI:MGI:",
      sapply(strsplit(identifiers, "MGI:"), tail, 1)
    )
  }
  # Map gene identifiers.
  suppressMessages({
    output <- AnnotationDbi::mapIds(osDB,
      keys = as.character(identifiers),
      column = columns(osDB)[colIDto],
      keytype = columns(osDB)[colIDfrom],
      multiVals = multiVals
    )
  })
  # Check if output is a list.
  if (is.list(output)) {
    # Replace NULL
    is_null <- sapply(output, is.null)
    output[is_null] <- NA
    output <- unlist(output)
  }
  # Check that all nodes (entrez) are mapped to gene symbols.
  not_mapped <- is.na(output)
  if (!quiet & sum(is.na(output)) != 0) {
    message(paste0(
      "Warning: Unable to map ", sum(not_mapped), " ", species, " ",
      from, "(s)", " to ", to, " identifiers!"
    ))
  }
  names(output) <- identifiers
  return(output)
}
