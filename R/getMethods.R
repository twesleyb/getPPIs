#' getMethods
#'
#' Given a HitPredict object (data.table), annotate the PPI-detection methods,
#' and refine PPI list by method type if provided.
#'
#' @param none
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
#' getMethods(hitpredict)
#' #
getMethods <- function(hitpredict, downloads = getwd(), keepdata = FALSE,
                       methods = "all", cutoff = 0.485) {
  # A method score >= 0.485 is considered to indicate high confidence (Villaveces et al.)
  # Download the Molecular Interactions Controlled Vocabulary.
  myfile <- file.path(downloads, "mi.owl")
  url <- "https://github.com/HUPO-PSI/psi-mi-CV/raw/master/psi-mi.obo"
  if (!file.exists(myfile)) {
    message("Downloading molecular ontology!")
    download.file(url, destfile = myfile)
    ontology <- ontologyIndex::get_ontology(myfile)
  } else {
    message("Using previously downloaded molecular ontology database!")
    ontology <- ontologyIndex::get_ontology(myfile)
  }
  # Remove downloaded files?
  if (keepdata == FALSE) {
    unlink(myfile)
  }
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
  # Filter by Confidence_score
  subdat <- hitpredict %>% filter(Confidence_score > cutoff)
  # Filter by methods.
  if (methods == "all") {
    return(hitpredict)
  } else {
    hitpredict <- hitpredict %>%
      dplyr::filter(Interaction_detection_methods %in% methods)
    return(hitpredict)
  }
}
