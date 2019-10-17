#' getProfile
#'
#' Partition the network with the Leiden algorithm.
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
#' getProfile()
getProfile <- function(g, partition_type, nsteps = 1, ...) {
  # Imports.
  suppressPackageStartupMessages({
    require(leiden)
  })
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
