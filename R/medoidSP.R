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
#' @import sarp.snowprofile
#'
#' @param profileList List of snowprofile objects
#' @param rescale_resample Do you want to uniformly rescale and resample the set of profiles prior to calculating the distance matrix?
#' @param retDistmat Do you want to *return* the pairwise distance matrix?
#' @param distmat If you have a precalculated distance matrix, provide it here to compute the medoid on it.
#' @param verbose print pairwise distance matrix? default FALSE
#' @param resamplingRate The resampling rate that is used for the whole set if `rescale_resample = TRUE`
#' @param ... arguments passed to [distanceSP]
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
#'   plot(grouplist, SortMethod = 'unsorted', labelOriginalIndices = TRUE)
#'
#'   ## calulate medoid profile
#'   idxMedoid <- medoidSP(grouplist)
#'   representativeProfile <- grouplist[[idxMedoid]]
#'   plot(representativeProfile, main = paste0("medoid (i.e., profile ", idxMedoid, ")"))
#'
#' }
#' @export

medoidSP <- function(profileList = NULL, rescale_resample = TRUE, retDistmat = FALSE, distmat = NULL, verbose = FALSE,
                     resamplingRate = 0.5, ...) {


  ## compute distance matrix from profileList:
  if (is.null(distmat)) {
    ## check input
    stopifnot(is.list(profileList))
    sapply(profileList, function(x) if (!is.snowprofile(x)) stop("At least one element in profileList is not a snowprofile"))

    ## rescale and resample
    if (rescale_resample) profileList <- reScaleSampleSPx(profileList, resamplingRate = resamplingRate)$set

    npros <- length(profileList)

    message(paste0("You are about to compute ", npros, "^2 = ", npros**2, " profile alignments. This will take roughly ",
               npros**2 * 0.3, " seconds, depending on the profile depths and resampling rate. Starting now..\n"))

    ## calculate pairwise dtw alignment between all profiles in profileList
    tstart <- Sys.time()
    pairDMat <- sapply(profileList, function(x)
      sapply(profileList, function(y)
        distanceSP(x, y, ...)
      )
    )
    message(paste0("Computed pairwise distance matrix. It actually took \n"))
    message(print(Sys.time()-tstart))

  } else {
    pairDMat <- distmat
  }

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
                 1, max)
  ## fill back the values twice, with transposing the matrix in between:
  pairDMat[iLT] <- symLT
  pairDMat <- t(pairDMat)
  pairDMat[iLT] <- symLT
  # compute sum of distances to all other profiles within group
  d <- rowSums(pairDMat)
  iMedoid <- which.min(d)

  if (retDistmat) return(list(iMedoid = iMedoid, distmat = pairDMat))
  else return(iMedoid)
}
