#' Average a group of snow profiles
#'
#' The functions [dbaSP] and [averageSP] implement Dynamic Time Warping Barycenter Averaging of snow profiles.
#' The convenient wrapper [averageSP] takes care of choosing several appropriate initial conditions and picking the optimal end result (by minimizing the mean squared error
#' between the average profile and the profile set). To pay appropriate attention to (thin) weak layers, weak layers need to be labeled in the profiles.
#' You can either do that manually before calling this routine to suit your personal needs, or you can provide specific properties (in `classifyPWLs`)
#' so that weak layers be labeled according to these properties by [sarp.snowprofile::labelPWL].
#' For more details, refer to the reference paper.
#'
#'
#' @describeIn averageSP convenient wrapper function
#' @param SPx SPx a [snowprofileSet] object. Note that the profile layers need to contain a column
#' called `$layerOfInterest` which classifies weak layers. While [averageSP] will label weak layers automatically if not done by the user beforehand, [dbaSP] won't do that but fail instead!;
#' consider thinking about how you want to label weak layers, see Description, `classifyPWLs` below, and the references.
#' Also note, that if you wish to average the *rescaled* profile set, do so manually before calling this function (see examples).
#' @param n the number of initial conditions that will be used to run [dbaSP]; see also [chooseICavg].
#' @param sm a [summary] of `SPx` metadata
#' @param progressbar should a progressbar be displayed (the larger n, the more meaningful the progressbar)
#' @param progressbar_pretext a character string to be prepended to the progressbar (mainly used by higher level cluster function)
#' @param classifyPWLs an argument list for a function call to [sarp.snowprofile::findPWL] which returns relevant PWLs for identifying initial conditions. **Importantly**, these arguments will also be used
#' to label weak layers in the profiles, if these labels do not yet exist in the layers objects as column `$layerOfInterest`.
#' Check out the documentation of `findPWL` to familiarize yourself with your manifold options!
#' @param classifyCRs an argument list for a function call to [sarp.snowprofile::findPWL] which returns relevant crusts for identifying initial conditions.
#' @param proportionPWL decimal number that specifies the proportion required to average an ensemble of grain types as weak layer type.
#' A value of 0.3, for example, means that layers will get averaged to a PWL type if 30% of the layers are of PWL type.
#' Meaningful range is between `[0.1, 0.5]`. Values larger than 0.5 get set to 0.5.
#' @param breakAtSim stop iterations when [simSP] between the last average profiles is beyond that value. Can range between `[0, 1]`. Default values differ between [dbaSP] and [averageSP].
#' @param breakAfter integer specifying how many values of simSP need to be above `breakAtSim` to stop iterating. Default values differ between [dbaSP] and [averageSP].
#' @param ... alignment configurations which are passed on to [dbaSP] and then further to [dtwSP]. Note, that you can't provide `rescale2refHS`, which is always set to FALSE. If you wish to rescale
#' the profiles, read the description of the `SPx` parameter and the examples.
#' @return A list of class `avgSP` that contains the fields
#'
#'   * `$avg`: the resulting average profile
#'   * `$set`: the corresponding resampled profiles of the group
#'   * `$call`: (only with `averageSP`) the function call
#'   * `$prelabeledPWLs`: (only with `averageSP`) boolean scalar whether PWLs (or any other layers of interest) were prelabeled before this routine (`TRUE`) or labeled by this routine (`FALSE`)
#'
#' @author fherla
#' @references Herla, F., Haegeli, P., and Mair, P.: Brief communication: A numerical tool for averaging large data sets of snow
#' stratigraphy profiles useful for avalanche forecasting, The Cryosphere Discuss., https://doi.org/10.5194/tc-2022-29, in review, 2022.
#'
#' @seealso [averageSPalongSeason]
#'
#' @examples
#' ## EXAMPLES OF averageSP
#' this_example_runs_about_10s <- TRUE
#' if (!this_example_runs_about_10s) {  # exclude from cran checks
#'
#' ## compute the average profile of the demo object 'SPgroup'
#' ## * by labeling SH/DH layers as weak layers,
#' ##   - choosing 3 initial conditions with an above average number of weak layers
#' ##   - in as many depth ranges as possible
#' ## * and neglecting crusts for initial conditions
#'
#'   avgList <- averageSP(SPgroup, n = 3,
#'                        classifyPWLs = list(pwl_gtype = c("SH", "DH")),
#'                        classifyCRs = NULL)
#'
#'   opar <- par(mfrow = c(1, 2))
#'   plot(avgList$avg, ymax = max(summary(avgList$set)$hs))
#'   plot(avgList$set, SortMethod = "unsorted", xticklabels = "originalIndices")
#'   par(opar)
#'
#'
#'   ## compute the average profile of the demo object 'SPgroup'
#' ## * by labeling SH/DH/FC/FCxr layers with an RTA threshold of 0.65 as weak layers,
#' ## * otherwise as above
#'
#'   SPx <- snowprofileSet(lapply(SPgroup, computeRTA))
#'   avgList <- averageSP(SPx, n = 3,
#'                        classifyPWLs = list(pwl_gtype = c("SH", "DH", "FC", "FCxr"),
#'                                            threshold_RTA = 0.65),
#'                        classifyCRs = NULL)
#'
#'   opar <- par(mfrow = c(1, 2))
#'   plot(avgList$avg, ymax = max(summary(avgList$set)$hs))
#'   plot(avgList$set, SortMethod = "unsorted", xticklabels = "originalIndices")
#'   par(opar)
#'
#' }
#'
#'
#' @export
averageSP <- function(SPx, n = 5, sm = summary(SPx),
                      progressbar = require("progress", quietly = TRUE, character.only = TRUE),
                      progressbar_pretext = NULL,
                      classifyPWLs = list(pwl_gtype = c("SH", "DH")),
                      classifyCRs = list(pwl_gtype = c("MFcr", "IF", "IFsc", "IFrc")),
                      proportionPWL = 0.5,
                      breakAtSim = 0.9, breakAfter = 2, verbose = FALSE,
                      ...) {

  if (!is.snowprofileSet(SPx)) stop("SPx must be a snowprofileSet")

  ## if the layer objects from SPx do not contain the column '$labeledPWL', label weak layers with the provided argument list
  prelabeledPWLs <- TRUE
  if (!all(sapply(SPx, function(sp) "layerOfInterest" %in% names(sp$layers)))) {
    message("Weak layers not pre-labeled. Automatic labeling based on arguments in classifyPWLs.")
    SPx <- snowprofileSet(lapply(SPx, function(sp) do.call("labelPWL", c(quote(sp), classifyPWLs))))
    prelabeledPWLs <- FALSE
  }

  ## compute several meaningful initial condition profiles:
  IC_ids <- chooseICavg(SPx, n = n, classifyPWLs = classifyPWLs, classifyCRs = classifyCRs, sm = sm)
  n <- length(IC_ids)
  sapply(seq(n), function(i) {
    if (!is.snowprofile(SPx[[IC_ids[i]]])) stop(paste0("Problem with choosing initial conditions! (idx: ", IC_ids[i], ")"))
  })

  ## initialize progressbar:
  if (progressbar) {
    pb <- progress::progress_bar$new(
      format = paste0(progressbar_pretext, " [:bar] :percent in :elapsed | eta: :eta"),
      total = n, clear = FALSE, width= 60)
  }

  ## compute average for different IC
  ## IC gets rescaled to median hs though!
  DBA <- lapply(seq(n), function(i) {
    tryCatch({
      dbaSP(SPx,
            scaleSnowHeight(SPx[[IC_ids[i]]], height = median(sm$hs))$queryScaled,  # IC profile rescaled to median hs
            sm = sm,
            proportionPWL = proportionPWL,
            breakAtSim = breakAtSim, breakAfter = breakAfter,
            verbose = verbose, ...)
    },
    error = function(err) {
      warning(paste0("Error in averaging of profiles:\n ", err))
      return(NA)
    }, finally = {if (progressbar) pb$tick()})
  })

  ## check whether there's a meaningful result:
  DBAmeaningful <- sapply(DBA, function(dba) ifelse(all(is.na(dba)), FALSE, TRUE))
  DBA <- DBA[DBAmeaningful]
  n <- length(DBA)
  if (n == 0) {
    if (median(sm$hs) < 2) stop("Can't find an average profile! Likely because your provided profiles are extremely shallow.")
    else stop("Can't find an average profile!")
  }

  ## compute MSE for each realization of average profile
  MSE <- sapply(seq(n), function(i) {
    DBA[[i]]$avg$mse
  })

  ## pick the average profile that has the smallest mse:
  ans <- DBA[[which.min(MSE)]]

  ## append other useful information:
  ans$call <- match.call()
  ans$prelabeledPWLs <- prelabeledPWLs



  return(ans)
}
