#' Calculate DTW alignment of two snow profiles
#'
#' This is the core function of the package and allows to match layers between pairs of snow profiles to align them. It
#' provides a variety of options, where the default values represent a good starting point to the alignment of most generic
#' profiles.
#'
#' @details
#' The individual steps of aligning snow profiles:
#'
#'  1. (optional) **Rescale** the profiles to the same height (cf., [scaleSnowHeight])
#'  2. **Resample** the profiles onto the same depth grid. 2 different approaches:
#'      - regular grid with a sampling rate that is provided by the user (recommended, cf., [resampleSP]).
#'      This approach requires to rescale the profiles.
#'      - irregular grid that includes all layer interfaces within the two profiles (i.e., set `resamplingRate = NA`) (cf., [resampleSPpairs])
#'  3. Compute a weighted **local cost matrix** from multiple layer characteristics (cf., [distMatSP])
#'  4. **Match the layers** of the profiles with a call to [dtw] (eponymous R package)
#'  5. Align the profiles by **warping** the query profile onto the reference profile (cf., [warpSP])
#'  6. (optional) If the function has been called with multiple different boundary conditions (global, top-down, or bottom-up alignments),
#'  the optimal alignment as determined by [simSP] will be returned.
#'
#'
#' @note
#' While DTW can be applied to sequences of different lengths, it is necessary that the two snow profiles have the same
#' number of layers for optimally warping one profile onto the other one and thereby aligning them. As long as the profiles
#' are rescaled and resampled this requirement is satisfied. This requirement is inconvenient only if the profiles
#' are not supposed to be rescaled (e.g. in a layer tracking application). In these cases, a workaround is to rescale the
#' profiles anyway, but increase the window size of the warping window proportionally to the snow height offset. In these cases
#' also consider including layer date information to reduce the likelihood of bad alignments.
#'
#' Furthermore, the alignment based on grain type information is currently only possible for specific grain types. These grain types
#' require a pre-defined distance or similarity, such as given by [grainSimilarity_align]. If your profile contains other grain types,
#' you are required to define your custom `grainSimilarity` matrix.
#'
#' @import dtw
#'
#' @param query The query snow profile to be warped
#' @param ref The reference snow profile to be warped against
#' @param open.end Is an open end alignment desired?
#' @param checkGlobalAlignment If `open.end = TRUE`, do you want to check whether a global alignment performs better (i.e.,
#' `open.end = FALSE`), and use the optimal one of the computed alignments?
#' @param keep.internals Append resampled and aligned snow profiles as well as internal parameters to the output object?
#' @param step.pattern The local slope constraint of the warping path, defaults to Sakoe-Chiba's symmetric pattern
#' described by a slope factor of P = 1, see [dtw::stepPattern]
#' @param resamplingRate Resampling rate for a regular depth grid; can also be a vector that provides the depth grid.
#' Set to `NA` to keep the (original, likely irregular) depth grid (see Details, bullet point 2.2).
#' @param rescale2refHS Rescale the query snow height to match the ref snow height?
#' @param bottom.up Compute an open.end alignment from the ground upwards?
#' @param top.down Compute an open.end alignment from the snow surface downwards?
#' @param ... Arguments passed to \code{\link{distMatSP}} and \code{\link{dtw}}, e.g.
#'
#'   * `dims`, `weights` (defaults specified in \code{\link{distMatSP}})
#'   * `ddateNorm`, numeric, normalize deposition date (default specified in \code{\link{distMatSP}})
#'   * `windowFunction`, default \code{\link{warpWindowSP}}
#'   * `window.size`, `ddate.window.size` (defaults specified in \code{\link{warpWindowSP}})
#'   * `gtype_distMat`, (default specified in \code{\link{distMatSP}}), cf. e.g. [grainSimilarity_align]
#'   * `prefLayerWeights`, weighting matrix for preferential layer matching, e.g. [layerWeightingMat]
#' @md
#'
#' @return
#' An alignment object of class 'dtwSP' is returned. This is essentially a list with various information about the alignment.
#' If `keep.internals = TRUE`, the resampled snow profiles 'query', 'reference' and 'queryWarped', as well as the
#' 'costMatrix' and 'directionMatrix' are elements of the returned object.
#'
#' @author fherla
#'
#' @seealso [plotSPalignment], [simSP]
#'
#' @examples
#'
#' dtwAlignment <- dtwSP(SPpairs$A_modeled, SPpairs$A_manual, open.end = FALSE)
#'
#' ## dtwSP object:
#' summary(dtwAlignment)
#' plot(dtwAlignment$queryWarped)
#' plotSPalignment(dtwAlignment = dtwAlignment)
#'
#' @export
dtwSP <- function(query, ref, open.end = TRUE, checkGlobalAlignment = TRUE, keep.internals = TRUE,
                  step.pattern = symmetricP1, resamplingRate = 0.5, rescale2refHS = TRUE,
                  bottom.up = TRUE, top.down = TRUE, ...) {

  ## --- assertion, setup etc ----
  if (!is.snowprofile(query) | !is.snowprofile(ref)) stop("query and ref need to be two snowprofile objects.")
  if (bottom.up == FALSE && top.down == FALSE) stop("Either bottom.up or top.down must be TRUE!")
  ## which direction to run dtw:
  dirs <- c("bottomUp", "topDown")[which(c(bottom.up, top.down))]
  if (checkGlobalAlignment && open.end) dirs <- c(dirs, "globalAlignment")
  ndirs <- length(dirs)
  if (ndirs > 1) keep.internals <- TRUE  # needed!

  ## --- rescaling ----
  if (rescale2refHS) {
    SC <- scaleSnowHeight(query, ref)
    query <- SC$queryScaled
    fac_trueQueryHeight <- SC$trueHeightFactor
  }

  ## --- resample profiles ----
  if (!is.na(resamplingRate)) {
    if (length(resamplingRate) == 1 & abs(query$hs - ref$hs) > resamplingRate) {
      ## one more check to rule out rounding errors and discrepancies in the hs naming:
      if (abs(max(query$layers$height) - max(ref$layers$height)) > resamplingRate) {
        stop("dtwSP: query and ref have different number of layers. Either set rescale2refHS to TRUE or provide a specific depth grid vector!")
      }
    }
    RES <- list(ref = resampleSP(ref, h = resamplingRate))
    RES$query <- resampleSP(query, h = RES$ref$layers$height)
  } else {
    RES <- resampleSPpairs(query, ref)
  }

  ## number of layers:
  nL <- nrow(RES$ref$layers)

  ## --- calculate dtw alignment ----
  DMat <- array(dim = c(nL, nL, ndirs))
  DMat[,,1] <- distMatSP(RES$query, RES$ref, ...)
  if (top.down) {
    ## mirror matrix at top-left/bottom-right diagonal:
    DMat[,, which(dirs == "topDown")] <- apply(apply(DMat[,,1], 1, rev), 1, rev)
  }
  if (checkGlobalAlignment) DMat[,, which(dirs == "globalAlignment")] <- DMat[,,1]
  A <- list()
  for (i in seq_along(dirs)) {
    if (dirs[i] %in% c("bottomUp", "topDown")) {
      ## Open End warps:
      ## OE warps are one-sided in the dtw package (i.e. only sub-sequences of the reference ae considered)
      ## in order to get fully symmetric OE warps, the roles of the reference and query are swapped and
      ## the winning pair is determined by the smallest normalized distance
      aTrue <- dtw(DMat[,,i], open.end = open.end, keep.internals = keep.internals, step.pattern = step.pattern, ...)
      ## for swapped roles transpose the distance matrix between the two profiles:
      aSwap <- dtw(t(DMat[,,i]), open.end = open.end, keep.internals = keep.internals, step.pattern = step.pattern, ...)
      ## choose winning pair:
      if (aTrue$normalizedDistance <= aSwap$normalizedDistance) {
        ## aTrue wins
        A[[i]] <- c(aTrue[c("costMatrix", "stepPattern", "N", "M", "openEnd", "openBegin", "windowFunction",
                                  "jmin", "distance", "normalizedDistance", "localCostMatrix")],
                          list(indexRef = aTrue$index2, index2 = aTrue$index2,  # indexRef is more clear, but keep index2 b/c it's implemented in functions
                               indexQuery = aTrue$index1, index1 = aTrue$index1,  # --=--
                               openEndType = "jmin", openEndType_verbose = "matched_subsequence_of_reference",
                               call = match.call(), direction = dirs[i])
        )
      } else {
        ## aSwap wins
        ## note, that we want to keep the original reference and query roles nonetheless
        ## i.e. - take true cost matrix (etc.)
        ##      - take distance from swap
        ##      - also take indices from swap, but swap them (i.e. mirrored at the diagonal); and jmin -> imin
        A[[i]] <- c(aTrue[c("costMatrix", "stepPattern", "N", "M", "openEnd", "openBegin", "windowFunction",
                                  "localCostMatrix")],
                          aSwap[c("distance", "normalizedDistance")],
                          list(indexRef = aSwap$index1, index2 = aSwap$index1,
                               indexQuery = aSwap$index2, index1 = aSwap$index2,
                               imin = aSwap$jmin, openEndType = "imin", openEndType_verbose = "matched_subsequence_of_query",
                               call = match.call(), direction = dirs[i])
        )
      }
    } else if (dirs[i] == "globalAlignment") {
      ## no OE warp / global alignment:
      A[[i]] <- dtw(DMat[,,i], keep.internals = keep.internals, open.end = FALSE, step.pattern = step.pattern, ...
      )[c("costMatrix", "stepPattern", "N", "M", "openEnd", "openBegin", "windowFunction",
          "jmin", "distance", "normalizedDistance", "localCostMatrix", "index1", "index2")]
      A[[i]] <- c(A[[i]], list(indexRef = A[[i]]$index2, indexQuery = A[[i]]$index1,
                                           openEndType = "jmin", open_EndType_verbose = "matched_full_sequences",
                                           call = match.call(), direction = dirs[i]))
    }  # END IF

    ## correct topDown alignment object (indices and matrices)
    if (A[[i]]$direction == "topDown") {
      A[[i]]$costMatrix <- apply(apply(A[[i]]$costMatrix, 1, rev), 1, rev)
      A[[i]]$localCostMatrix <- apply(apply(A[[i]]$localCostMatrix, 1, rev), 1, rev)
      A[[i]]$index1 <- A[[i]]$N + 1 - A[[i]]$index1
      A[[i]]$index2 <- A[[i]]$M + 1 - A[[i]]$index2
      A[[i]][A[[i]]$openEndType] <-  A[[i]]$N + 1 - A[[i]][A[[i]]$openEndType][[1]]
    }
    ## append profiles to alignment object
    if (keep.internals) {
      A[[i]]$query <- RES$query
      A[[i]]$reference <- RES$ref
      A[[i]] <- warpSP(A[[i]])
      if (rescale2refHS) A[[i]]$fac_trueQueryHeight <- fac_trueQueryHeight
    }
  }  # END LOOP

  ## choose whether bottomUp / TopDown / globalAlignment performed better (based on external function call simSP):
  if (ndirs > 1) {
    gtDM <- sim2dist(grainSimilarity_evaluate(FALSE))
    sim <- rep(NA, times = ndirs)
    for (i in seq_along(dirs)) {
      A[[i]]["sim"] <- sim[i] <-  simSP(RES$ref, A[[i]]$queryWarped, gtype_distMat = gtDM)
    }
    win <- which.max(sim)
  } else {
    win <- 1
  }

  ## modify local cost matrix to resemble step pattern constraints:
  A[[win]]$localCostMatrix[is.na(A[[win]]$costMatrix)] <- NA

  ## inherits new class due to changes made to object, but keeps dtw class to use existing functionality (esp. dtwPlotDensity())
  class(A[[win]]) <- append("dtwSP", class(A[[win]]))



  return(A[[win]])

}
