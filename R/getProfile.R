#' getProfile
#'
#' Partition the network with the Leiden algorithm.
#'
#' @param g - igraph object
#'
#' @param partition_type - Leiden algorithm method to be used.
#'
#' @param nsteps - number of resolutions to be analyzed. If 1 then,
#' the native resolution (r=0) will be used.
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
#' getProfile()
getProfile <- function(g, partition_type, nsteps = 1, ...) {
  # Imports.
  suppressPackageStartupMessages({
    require(leiden)
  })
  # Initialize progress bar.
  if (i == 1) {
    message(paste("Examining optimal partitioning of the graph at", nsteps, "resolutions."))
    pb <- txtProgressBar(min = 0, max = nsteps, initial = 0)
  }
  setTxtProgressBar(pb, i)
  # Loop through resolutions and perform leiden algorithm clustering using specified method.
  resolution <- seq(0, 1, length.out = nsteps)
  profile <- list()
  for (i in 1:nsteps) {
    partition <- leiden(g,
      partition_type,
      resolution_parameter = resolution[i],
      seed = NULL,
      n_iterations = -1
    )
    profile[[i]] <- partition
  }
  return(profile)
}
