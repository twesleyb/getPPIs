#' installDBs
#'
#' Installation of organism specific mapping databases from BioConductor.
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
#' @import BiocManager
#'
#' @export
#'
#' @examples
#' installDBs()
installDBs <- function() {
  # Fixme: should check if package already exists!
  annotationDBs <- mappingDBs()
  # Loop to install packages.
  packages <- unlist(sapply(annotationDBs, "[", 2))
  for (package in packages) {
    if (requireNamespace(package, quietly = TRUE)) {
      message(paste(package, "is already installed!"))
    } else {
      BiocManager::install(package, update = FALSE, ask = FALSE)
    }
  }
}
