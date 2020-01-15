#' getHomologs
#'
#' Downloading the NCBI holoogy database. Mapping proteins to species of
#' interest.
#'
#' @param hitpredict (data.table) - HitPredict data.
#'
#' @param taxid (integer) - taxonomic identifier for species of interest
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
#' getHomologs(entrez, taxid = 10090)
getHomologs <- function(entrez, taxid, verbose = FALSE) {
  suppressMessages({
    require(getPPIs)
  })
  # Download and load NCBI homology gene data.
  downloads <- getwd()
  url <- "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"
  destfile <- file.path(downloads, "homologene.data")
  if (verbose) {
    message("Downloading NCBI HomoloGene data...")
  }
  download.file(url, destfile, quiet = !verbose)
  homologene <- data.table::fread(destfile, header = FALSE)
  # Remove raw data.
  unlink(destfile)
  # Fix column names.
  # Gene ID is a genes organism specific Entrez ID.
  # HID is the genes homology id.
  homologene <- homologene %>% rename_all(list(~ c(
    "HID", "TaxonomyID", "GeneID",
    "GeneSymbol", "ProteinGI", "ProteinAccession"
  )))
  # Use HomoloGene to create homology map.
  homology_map <- as.list(homologene$HID)
  names(homology_map) <- homologene$GeneID
  # Map Entrez to homology ID.
  hid <- homology_map[entrez]
  hid[seq(1:length(hid))[unlist(lapply(hid, is.null))]] <- NA
  # Subset homologene data, keep all genes from your species of interest.
  homologene <- homologene %>% filter(TaxonomyID == taxid)
  # Taxonomy info.
  annotationDBs <- mappingDBs()
  osDB <- unlist(annotationDBs[sapply(annotationDBs, "[", 1) == taxid])
  names(osDB) <- sapply(strsplit(names(osDB), "\\."), "[", 2)
  organism <- osDB["alias"]
  # Get entrez associated with these HIDs.
  osEntrez <- homologene$GeneID[match(hid, homologene$HID)]
  # Status report.
  is_missing <- is.na(osEntrez)
  n <- length(is_missing)
  if (verbose) {
    message(paste0(
      round(100 * sum(!is_missing) / n, 2),
      "% of supplied entrez IDs were mapped to homologous genes in ",
      organism, "."
    ))
  }
  # osEntrez <- osEntrez[!is_missing]
  # Return data with genes mapped to gene specific entrez.
  return(osEntrez)
}
