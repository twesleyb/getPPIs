#' getInteractionMethods
#'
#' Given a HitPredict object (data.table), annotate the PPI-detection methods.
#' Note: A method score >= 0.485 is considered to indicate high confidence
#' according to Villaveces et al.
#'
#' @param hitpredict (data.table)
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords annotation database org.db entrez
#'
#' @export
#'
#' @examples
#' getInteractionMethods(hitpredict)
getInteractionMethods <- function(hitpredict,quiet=TRUE) {
  # Download the molecular ontology.
  downloads <- getwd()
  myfile <- file.path(downloads, "mi.owl")
  url <- "https://github.com/HUPO-PSI/psi-mi-CV/raw/master/psi-mi.obo"
  if (!quiet) { message("Downloading molecular interaction ontology!") }
  download.file(url, destfile = myfile,quiet=quiet)
  ontology <- ontologyIndex::get_ontology(myfile)
  unlink(myfile)
  # Function to extact method names.
  getMethod <- function(x) {
    y <- gsub("psi-mi:", "", unlist(strsplit(x, "\\|")))
    z <- ontology$name[y]
    return(z)
  }
  # Get method names.
  mi <- hitpredict$Interaction_detection_methods
  meth <- lapply(mi, getMethod)
  namen <- unlist(lapply(meth, function(x) paste(x, collapse = " | ")))
  # Add this to hitpredict.
  hitpredict <- hitpredict %>% dplyr::mutate(Methods = namen)
  # Return data.
  return(hitpredict)
}
