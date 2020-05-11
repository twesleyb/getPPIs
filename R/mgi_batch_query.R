#' mgi_batch_query()
#'
#' map uniprot IDs to entrez gene IDs using the MGI batch query tool.
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
#' @export
#'
mgi_batch_query <- function(ids, from = "Uniprot", to = "Entrez",
                            download = TRUE, quiet = TRUE) {

  # Map gene identifiers using MGI mapping data.
  # FIXME: Currently only supports uniprot to entrez mapping.

  # Imports
  suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
  })

  # Download the MGI mapping data.
  # Grab the gene identifier columns.
  # Separate rows with multiple gene identifiers.
  # Remove rows with no matching Uniprot ID.
  # Collect MGI and Uniprot IDs in a data.table.

  # "http://www.informatics.jax.org/downloads/reports/MRK_Sequence.rpt"
  url <- "https://bit.ly/3et6O4M"

  # All gene identifiers.
  IDs <- c(
    "MGI Marker Accession ID", "GenBank IDs",
    "Ensembl transcript IDs",
    "TrEMBL IDs", "RefSeq protein IDs", "Marker Symbol",
    "RefSeq transcript IDs", "UniProt IDs", "Ensembl protein IDs",
    "UniGene IDs"
  )

  # Get column cooresponding to input gene identifier type (from).
  idy <- grep(tolower(from), tolower(IDs))
  geneID <- IDs[idy]
  if (length(geneID) > 1) {
    stop("Multiple matching gene identifiers.")
  }

  # Load MGI gene mapping data.
  if (!download) {
    data(MGI_Gene_Map)
  } else {
    # Build a mappint database.
    message("Creating mapping database...")
    # Download the data.
    MRK_Sequence <- fread(url, showProgress = !quiet)
    # Get MGI and gene ID columns.
    cols <- c("MGI Marker Accession ID", geneID)
    dtA <- MRK_Sequence %>%
      dplyr::select(all_of(cols))
    # Simplify column names.
    colNames <- c(
      "MGI", "GenBank", "Ensembl_transcript", "TrEMBL",
      "RefSeq_protein", "Symbol", "RefSeq_transcript",
      "UniProt", "Ensembl_protein", "UniGene"
    )
    colnames(dtA) <- c("MGI", colNames[idy])
    colName <- colNames[idy]
    # Dynamically execute command to separate rows.
    cmd <- paste0(
      "tidyr::separate_rows(dtA,",
      colName, ",sep='\\\\|')"
    )
    dtA <- eval(parse(text = cmd))
    # Dynamically execute command to remove empty rows.
    cmd <- paste0("filter(dtA,", colName, "!='')")
    dtA <- eval(parse(text = cmd))
    # Coerce result to a data.table.
    dtA <- as.data.table(dtA)
    # Download MGI to Entrez mapping data.
    # Separate rows with multiple secondary MGI IDs.
    # Collect Entrez and MGI ids in a data.table.
    # http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt
    url <- "https://bit.ly/3aeO1Xi"
    MGI_EntrezGene <- fread(url, showProgress = !quiet)
    dtB <- MGI_EntrezGene %>% dplyr::select(V1, V8, V9)
    colnames(dtB) <- c("MGI_A", "MGI_B", "Entrez")
    dtB <- tidyr::separate_rows(dtB, "MGI_B", sep = "\\|")
    dtB <- melt(dtB,
      id.vars = "Entrez",
      measure.vars = c("MGI_A", "MGI_B"), value.name = "MGI"
    )
    dtB <- dtB %>%
      dplyr::select(Entrez, MGI) %>%
      na.omit()
    # Merge the two tables.
    dtC <- inner_join(dtA, dtB, by = "MGI") %>% as.data.table()
    MGI_Gene_Map <- dtC
  } # Ends chunk that builds a map.

  # Map IDs.
  idx <- match(ids, MGI_Gene_Map$UniProt)
  entrez <- MGI_Gene_Map$Entrez[idx]
  n_missing <- sum(is.na(entrez))
  if (n_missing > 0) {
    warn(paste(
      "Unable to map", n_missing,
      "gene identifiers to entrez."
    ))
  }
  names(entrez) <- ids

  # Return entrez ids.
  return(entrez)
}
