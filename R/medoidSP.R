#' Find the medoid snow profile among a group of profiles
#'
#' Find the medoid snowprofile among a group of profiles, based on their pairwise dissimilarity. Either provide a list
#' of `snowprofile` objects, or a precomputed distance matrix. \cr
#' If you provide a list of profiles the profiles can optionally be rescaled and resampled before the distance matrix
#' for the medoid calculation is computed. When computing the distance matrix this routine calls [distanceSP] for
#' *every possible pair* of profiles among the group. During that call the profile pair is aligned by [dtwSP]
#'  and the aligned pair is evaluated by [simSP].
#' Note that the number of possible profile pairs grows exponentially with the number of profiles in the group (i.e.,
#' O(n^2) calls, where n is the number of profiles in the group).
#'
#' Note that the pairwise distance matrix is modified within the function call to represent a symmetric distance matrix.
#' That is,, however, not originally the case, since `dtwSP(A, B) != dtwSP(B, A)`. The matrix is therefore made symmetric by
#' setting the similarity between the profiles A and B to `max({dtwSP(A, B), dtwSP(B, A)})`.
#'
#' @import sarp.snowprofile
#'
#' @param profileList List of snowprofile objects
#' @param rescale_resample Do you want to uniformly rescale and resample the set of profiles prior to calculating the distance matrix?
#' @param retDistmat Do you want to *return* the pairwise distance matrix?
#' @param distmat If you have a precalculated distance matrix, provide it here to compute the medoid on it.
#' @param verbose print pairwise distance matrix? default FALSE
#' @param resamplingRate The resampling rate that is used for the whole set if `rescale_resample = TRUE`
#' @param progressbar Do you want to print a progress bar with recommended package "progress"?
#' @param ... arguments passed to [distanceSP] and then further to [dtwSP]
#'
#' @return If `retDistmat = FALSE` return the (named) index of the medoid snow profile, otherwise return a list with the elements
#' `iMedoid` and `distmat`.
#'
#' @author fherla
#' @seealso [reScaleSampleSPx]
#'
#' @examples
#' this_example_runs_about_5s <- TRUE
#' if (!this_example_runs_about_5s) {  # exclude from cran checks
#'
#'   ## take a list of profiles
#'   grouplist <- SPgroup[1:4]
#'   plot(grouplist, SortMethod = 'unsorted', xticklabels = "originalIndices")
#'
#'   ## calulate medoid profile
#'   idxMedoid <- medoidSP(grouplist)
#'   representativeProfile <- grouplist[[idxMedoid]]
#'   plot(representativeProfile, main = paste0("medoid (i.e., profile ", idxMedoid, ")"))
#'
#' }
#' @export

medoidSP <- function(profileList = NULL, rescale_resample = TRUE, retDistmat = FALSE, distmat = NULL, verbose = FALSE,
                     resamplingRate = 0.5, progressbar = require("progress", quietly = TRUE, character.only = TRUE), ...) {


  ## compute distance matrix from profileList:
  if (is.null(distmat)) {
    ## check input
    stopifnot(is.list(profileList))
    sapply(profileList, function(x) if (!is.snowprofile(x)) stop("At least one element in profileList is not a snowprofile"))

    ## rescale and resample
    if (rescale_resample) {
      profileList <- reScaleSampleSPx(profileList, resamplingRate = resamplingRate)$set
    }

    npros <- length(profileList)

    if (npros > 10)
      message(paste0("You are about to compute ", npros, "**2 = ", npros**2, " profile alignments. Be patient.."))


    ## initialize progressbar:
    if (progressbar) {
      pb <- progress::progress_bar$new(
        format = " [:bar] :percent in :elapsed | eta: :eta",
        total = npros**2, clear = FALSE, width= 60)
    }

    ## the function won't fail upon errors, but ensure that these are printed as warnings right away upon occuring:
    op <- options("warn")
    on.exit(options(op))
    options(warn=1)

    ## initialize subfunction to calculate pairwise dtw alignment between all profiles in profileList
    ## with progressbar
    if (progressbar) {
      start_computation <- function(profileList, ...) {
        pairDMat <- sapply(profileList, function(x)
          sapply(profileList, function(y)
            tryCatch({
              pb$tick()
              distanceSP(x, y, ...)
            },
            error = function(err) {
              warning(paste0("Error in alignment of ", x$station_id, ", ", y$station_id, ":
                             ", err))
              return(NA)
            })
          )
        )
        return(pairDMat)
      }
    } else {
      ## and without progressbar
      start_computation <- function(profileList, ...) {
        pairDMat <- sapply(profileList, function(x)
          sapply(profileList, function(y)
            tryCatch({distanceSP(x, y, ...)},
                     error = function(err) {
                       warning(paste0("Error in alignment of  ", x$station_id, ", ", y$station_id, ":"))
                       warning(err)
                       return(NA)
                     })
          )
        )
        return(pairDMat)
      }
    }

    ## Do computations:
    pairDMat <- start_computation(profileList, ...)


  } else {
    pairDMat <- distmat
  }

  diag(pairDMat) <- 0  # ensure that auto-alignments are no partial alignments
  rownames(pairDMat) <- seq(nrow(pairDMat))
  colnames(pairDMat) <- seq(nrow(pairDMat))
  if (verbose) print(pairDMat)

  ## medoid element:
  ## (1) pairDMat is yet an asymmetric matrix [b\c dtw(a, b) != dtw(b, a)]
  ##     hence, make it symmetric by filling in max(element_in_lower_triangle, element_in_upper_triangle)
  ## (2) compute row sums and search for the (intra-group) minimal total distance, i.e. medoid
  iLT <- lower.tri(pairDMat)
  symLT <- apply(matrix(data = c(pairDMat[iLT],
                                 t(pairDMat)[iLT]),
                        ncol = 2),
                 1, max, na.rm = TRUE)
  ## fill back the values twice, with transposing the matrix in between:
  pairDMat[iLT] <- symLT
  pairDMat <- t(pairDMat)
  pairDMat[iLT] <- symLT
  # compute sum of distances to all other profiles within group
  d <- rowSums(pairDMat, na.rm = TRUE)
  iMedoid <- which.min(d)

  if (retDistmat) return(list(iMedoid = iMedoid, distmat = pairDMat))
  else return(iMedoid)
}
