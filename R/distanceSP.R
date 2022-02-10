#' Wrapper for dtwSP and simSP
#'
#' Calculate the distance between two snowprofile objects by
#'
#'  1. Matching their layers and aligning them (i.e., warp one profile onto the other one)
#'  2. Assessing the similarity of the aligned profiles based on avalanche hazard relevant characteristics
#'  3. Convert the similarity score into a distance value between `[0, 1]`
#'
#' @param query The query snowprofile object (will be warped onto `ref`)
#' @param ref The reference snowprofile object (will *not* be warped)
#' @param ... passed on to [dtwSP]
#'
#' @author fherla
#'
#' @seealso [dtwSP], [simSP], [medoidSP]
#'
#' @details This procedure is useful for clustering and aggregating tasks, given a set of multiple profiles.
#'
#' @export
distanceSP <- function(query, ref, ...) {


  ## first: DTW alignment of snow profiles
  ## second: compute distance via similarity measure for snow profiles
  distance <- tryCatch({
    (1 - dtwSP(query, ref, ...)$sim)
  }, error = function(err) {
    warning(err)
    1  # in case of alignment error, set distance to 1
    }
  )

  return(distance)
}
