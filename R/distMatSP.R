#' Calculate a multidimensional distance matrix between two profiles
#'
#' This routine calculates a distance matrix for two given profiles (`query` and `ref`). Analogously to other DTW
#' routines, the query is arranged along the matrix rows, the ref along the columns. Every cell of the matrix represents
#' the distance between the corresponding profile layers. The distance is calculated based on the specified layer properties
#' (e.g., `hardness`, `gtype`, `ddate`). The routine calls subroutines to calculate the distance for each property and
#' combines the normalized distances by weighted averaging.
#'
#' @param query The query snowprofile object
#' @param ref The ref snowprofile object
#' @param dims Character vector containing the layer properties to calculate the distance over. Currently implemented
#' are the properties `hardness`, `gtype`, `ddate`.
#' @param weights Numeric vector of the same length as `dims` specifying the averaging weights to each element of dims.
#' @param gtype_distMat A symmetric **distance** scoring matrix provided as data.frame that stores information about
#' the distances between the encountered grain types of the provided profiles. Default is the corresponding distance
#' matrix of [grainSimilarity_align], cf. [sim2dist].
#' @param prefLayerWeights A matrix similar to `gtype_distMat`, but storing weights for preferential layer matching,
#' e.g. defaults to [layerWeightingMat]; the higher the values for a given grain type pair, the more the algorithm will try to
#' match those layers above others. To turn weighting scheme off, set `prefLayerWeights = NA`
#' @param ddateNorm Normalize the deposition date distance by `ddateNorm` number of days. Numeric, default 5.
#' @param windowFunction a window function analogous to [warpWindowSP] (*Note!* designed for a quadratic distance matrix,
#' i.e. two profiles with identical numbers of layers. To resample the profiles see [resampleSPpairs] and examples below.
#' Other compatible window functions can be found in [dtw::dtwWindowingFunctions].)
#' @param ... arguments to the window function, e.g. `window.size`, `ddate.window.size`, ...
#'
#' @return A distance matrix of dimension (n x m), where n, m are the number of layers in the query and ref, respectively.
#'
#' @author fherla
#'
#' @seealso [resampleSPpairs]
#'
#' @note For package developers: dot inputs to the function (i.e., `...`) also necessary to keep [dtwSP] highly flexible
#' and customizable. Dot inputs may contain arguments that remain unused in this function.
#'
#' @examples
#'
#' ## call function with two snow profiles of unequal lengths, without using a window function:
#' dMat_noWindow <- distMatSP(SPpairs$A_modeled, SPpairs$A_manual, windowFunction = NA)
#' graphics::image(dMat_noWindow, main = "Default distance matrix without a warping window")
#'
#'
#' ## compute distance based on grain type alone,
#' ## and additionally disable preferential layer matching:
#' dMat <- distMatSP(SPpairs$A_modeled, SPpairs$A_manual, windowFunction = NA,
#'                   dims = "gtype", weights = 1, prefLayerWeights = NA)
#' graphics::image(dMat,
#'                 main = "... only based grain type, and without preferential layer matching")
#'
#'
#' ## to use a warping window, the profiles need to have equal numbers of layers:
#' ## i.e., resample them first, then compute dMat with a warping window
#' plist <- resampleSPpairs(SPpairs$A_modeled, SPpairs$A_manual)
#' dMat <- distMatSP(plist$query, plist$ref)
#' graphics::image(dMat, main = "Default with warping window")
#'
#' @export

## TODO: apply window function before calculating the distances of all elements to save a few resources

distMatSP <- function(query, ref, dims = c("hardness", "gtype"), weights = c(0.2, 0.8),
                                 gtype_distMat = sim2dist(grainSimilarity_align(FALSE)),
                                 prefLayerWeights = layerWeightingMat(FALSE),
                                 ddateNorm = 5, windowFunction = warpWindowSP, ...) {

  ## --- Assertions and Initializations ----
  if (!is.snowprofile(query) | !is.snowprofile(ref)) stop("Both query and ref need to be snowprofile objects!")
  nQue <- nrow(query$layers)
  nRef <- nrow(ref$layers)
  if (!length(dims) == length(weights)) stop("Weights and dims need to have equal length!")
  if (!all(dims %in% names(ref$layers)) & !all(dims %in% names(query$layers))) stop("Not all dimensions you specified are available in both profiles.")

  ## handle to pass to subfunctions for error handling:
  handle <- paste("Query @", query$station, query$datetimeUTC, "&",
                  "Ref @", ref$station, ref$datetimeUTC)

  ## initialize layer weighting matrix for preferential alignment
  lwmat <- array(0, c(nQue, nRef))

  ## initialize matrix with query along rows and reference along columns,
  ## third dimension is as long as dims
  distArr <- array(data = NA, dim = c(nQue, nRef, length(dims)))
  sQue <- seq(nQue)
  sRef <- seq(nRef)
  # Index matrix of dim (nQue x nRef)x3, that holds the index position of every field in distArr:
  iMat <- matrix(c(rep(sQue, each = nRef), rep(sRef, times = nQue), rep(NA, times = nQue*nRef)), ncol = 3)

  ## --- Calculate distances and store them in multi-dim array ----
  tofill <- 1  # which dimension needs to be filled next?
  if ("hardness" %in% dims) {
    iMat[,3] <- tofill
    distArr[iMat] <-
      weights[which(dims == "hardness")] *
      hardnessDistance(query$layers$hardness[iMat[,1]],
                       ref$layers$hardness[iMat[,2]],
                       normalize = TRUE,
                       absDist = TRUE)

    if (any(is.na(distArr[,,tofill]))) stop("NAs encountered in hardness distance")
    tofill <- tofill + 1
  }
  if ("gtype" %in% dims) {
    iMat[,3] <- tofill
    distArr[iMat] <-
      weights[which(dims == "gtype")] *
      extractFromScoringMatrix(ScoringFrame = gtype_distMat,
                               grainType1 = as.character(query$layers$gtype)[iMat[,1]],
                               grainType2 = as.character(ref$layers$gtype)[iMat[,2]],
                               profile_handle = handle)
    tofill <- tofill + 1

    ## fill layer weighting matrix
    if (!all(is.na(prefLayerWeights))) {
      lwmat[iMat[, 1:2]] <-
        extractFromScoringMatrix(ScoringFrame = prefLayerWeights,
                                 grainType1 = as.character(query$layers$gtype)[iMat[,1]],
                                 grainType2 = as.character(ref$layers$gtype)[iMat[,2]])
    }
  }
  if ("ddate" %in% dims) {
    iMat[,3] <- tofill
    distArr[iMat] <-
      weights[which(dims == "ddate")] *
      ddateDistance(query$layers$ddate[iMat[,1]],
                    ref$layers$ddate[iMat[,2]],
                    normalizeBy = ddateNorm)
    tofill <- tofill + 1
  }

  ## --- Reduce multi-dim array to final distance matrix ----
  ## Sum over array dimensions and multiply by layer weighting scheme to get one (nQue x nRef) distance matrix:
  if (any(is.na(distArr))) warning("Ignoring encountered NAs in distance calculation!")
  maxLW <- suppressWarnings(max(c(max(prefLayerWeights, na.rm = TRUE), 1)))
  distMat <- rowSums(distArr, dims = 2, na.rm = TRUE) + lwmat

  ## Apply window function to distMat:
  ## windowFunction needs to return TRUE for the warping window, FALSE else
  ## therefore, set all FALSE returns of the window function to NA in the distMat:
  if (is.function(windowFunction)) {
    if (nQue != nRef) stop("Warping window is only implemented for profiles with identical numbers of layers, resample them first! ..see ?distMatSP")
    profile.size <- nrow(ref$layers)
    profile.height <- ref$layers$height[profile.size]
    iheight <- matrix(query$layers$height, byrow = FALSE, nrow = profile.size, ncol = profile.size)
    jheight <- matrix(ref$layers$height, byrow = TRUE, nrow = profile.size, ncol = profile.size)
    if ("ddate" %in% names(query$layers) & "ddate" %in% names(ref$layers)) {
      iddate <- matrix(query$layers$ddate, byrow = FALSE, nrow = profile.size, ncol = profile.size)
      jddate <- matrix(ref$layers$ddate, byrow = TRUE, nrow = profile.size, ncol = profile.size)
    } else {
      iddate <- NA
      jddate <- NA
    }
    distMat[!windowFunction(row(distMat), col(distMat),
                            iheight = iheight, jheight = jheight,
                            iddate = iddate, jddate = jddate,
                            profile.size = profile.size, profile.height = profile.height, ...)] <- NA
  }

  ## final assertion (during code writing, can be eliminated if stable)
  if (!all(dim(distMat) == c(nrow(query$layers), nrow(ref$layers)))) stop("Calculated distance matrix seems to have wrong dimensions")

  return(distMat)
}
