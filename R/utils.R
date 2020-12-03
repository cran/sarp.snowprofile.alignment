##
## This file 'utils_alignment.R' contains useful helper routines for sarp.snowprofile.alignment,
## some of which are NOT exported, i.e. only available internally, but not available after
## attaching the package.
##
## Separate each function with a section.
##
## fherla
##



## --- Difference in Hand Hardness ----
#' Difference in Hand Hardness
#'
#' Calculate the difference (i.e. distance) in hand hardness
#'
#' @param hardness1 character or numeric hand hardness value (1D array)
#' @param hardness2 character or numeric hand hardness value (1D array)
#' @param normalize Should result be normalized? boolean, default False.
#' @param absDist Interested in absolute distance? default True.
#' @return numeric Hand Hardness Distance
#' @author fherla
hardnessDistance <- function(hardness1, hardness2, normalize=FALSE, absDist=TRUE) {

  ## Convert character HHI to numeric HHI if necessary
  if (!any(is.numeric(hardness1))) hn1<-unname(sapply(hardness1, char2numHHI)) else hn1<-hardness1
  if (!any(is.numeric(hardness2))) hn2<-unname(sapply(hardness2, char2numHHI)) else hn2<-hardness2

  ## Calculate the difference between vector entries
  hdist <- hn1 - hn2
  if (absDist) hdist <- abs(hdist)
  if (normalize) hdist <- hdist / 5

  return(hdist)
}

## --- Convert similarity to distance ----
#' Convert 'similarity' matrix to 'distance' matrix
#'
#' Convert a 'similarity' matrix to 'distance' matrix. *Note* that the similarity must be normalized (i.e. within \[0, 1\])
#' @param SimMat similarity matrix of type data.frame with ranges \[0, 1\]
#' @return copy of input data.frame with similarities inverted to distances (i.e. dist = 1 - sim)
#' @author fherla
#' @examples
#'
#' ## the 'swissSimilarityMatrix' as similarity and as distance
#' graphics::image(as.matrix(swissSimilarityMatrix))
#' graphics::image(as.matrix(sim2dist(swissSimilarityMatrix)))
#'
#' @export
sim2dist <- function(SimMat) {
  DistMat <- data.frame(SimMat)
  convert <- function(x) ifelse(is.na(x), NA, 1-x)
  DistMat[] <- lapply(DistMat, function(x) convert(x))
  return(DistMat)
}


## --- Extract from Scoring matrix ----
#' Extract from Scoring matrix
#'
#' Vectorized function to efficiently extract elements from scoring matrix of type data.frame
#' @param ScoringFrame Scoring matrix of type data.frame (needs to be of symmetric, matrix like format)
#' @param grainType1 character vector (yes, vector!) of grain type contained in ScoringFrame
#' @param grainType2 same as `grainType1`
#' @param profile_handle character or numeric handle that links a potential warning message to the set of grain types,
#' if an unknown grain type is encountered (must be of length = 1)
#' @return numeric vector of length `grainType1` with the elements of `ScoringFrame`
#' that are defined by `grainType1` and `grainType2`
#' @author fherla
#'
extractFromScoringMatrix <- function(ScoringFrame, grainType1, grainType2,
                                     profile_handle = NULL) {
  if (!all(ScoringFrame[upper.tri(ScoringFrame)] ==
           t(ScoringFrame)[upper.tri(ScoringFrame)])) {
    stop("need symmetric matrix-like dataframe!
         upper and lower triangle are not equal.")
  }
  if (any(is.na(grainType1)) | any(is.na(grainType2))) {
    warning(paste0("Missing grain types in profile ", profile_handle,
                   ", assigning default value as specified in ScoringFrame"))
  }
  gT1 <- ifelse(is.na(grainType1), "na", grainType1)
  gT2 <- ifelse(is.na(grainType2), "na", grainType2)

  d <- ScoringFrame[cbind(gT1, gT2)]

  return(d)
}


## --- Deposition Date Distance ----
#' Deposition Date Distance
#'
#' Calculate the distance (i.e. dissimilarity) between two deposition dates
#' @param ddate1 1D array of POSIX dates
#' @param ddate2 same format and length as ddate1
#' @param normalizeBy Numeric scalar to be used for normalization, i.e. the number of days, that defines the distance value of 1
#' @param clipWindow Should differences larger than 'normalizeBy' number of days be set to distance 'Infinity'? default FALSE.
#' @param na.dist replace NA values with that distance
#'
#' @return An array of length(ddate1) containing the distances according to the configurations.
#' @author fherla
#' @examples
#' ## create ddate arrays..
#' ddate <- as.POSIXct("2019/04/20 12:00", tz = "UTC")
#' ddate1 <- rep(ddate, 5)
#' ddate2 <- as.POSIXct(c("2019/04/12 08:00", "2019/04/16 10:00", "2019/04/20 12:00",
#'                        "2019/04/21 16:00", "2019/04/22 20:00"), tz = "UTC")
#'
#' ## .. and calculate distance:
#' ddateDistance(ddate1, ddate2, normalizeBy = 5)
#' @export
ddateDistance <- function(ddate1, ddate2, normalizeBy = 5, clipWindow = FALSE, na.dist = 0.5) {

  # convert from POSIX format to julian date number:
  jd1 <- as.numeric(julian(ddate1))
  jd2 <- as.numeric(julian(ddate2))
  # jd1 <- ifelse(is.na(ddate1), NA, as.numeric(julian(ddate1)))
  # jd2 <- ifelse(is.na(ddate1), NA, as.numeric(julian(ddate2)))
  ## Calculate the difference between vector entries
  jdist <- abs(jd1 - jd2)
  jdist <- jdist / normalizeBy
  jdist[is.na(jdist)] <- na.dist
  if (clipWindow) jdist[jdist > 1] <- Inf

  return(jdist)
}


## --- Height Scaling of Profiles ----
#' Scale total height of a snow profile
#'
#' Scale the snow height of a snow profile either (1) based on another profile, or (2) based on a provided (predetermined) snow height.
#' This function can therefore be used to scale two snow profiles to an identical snow height by scaling the height vector of the (query) profile
#' against the height vector of the (reference) profile.
#' @param query the query snow profile (whose height vector will be scaled)
#' @param ref the reference snow profile (whose total snow height will be used as the reference height for the scaling)
#' @param height an optional reference height that can be given instead of the query profile
#' @return query profile with scaled height vector
#' @author fherla
#' @export
scaleSnowHeight <- function(query, ref = NA, height = NA) {

  if (!is.snowprofile(query)) stop("query needs to be a snowprofile object")
  if (!all(is.na(ref))) {
    if(!is.snowprofile(ref)) stop("ref needs to be a snowprofile object")
    fac <- max(ref$layers$height)/max(query$layers$height)
  } else if (!is.na(height)) {
    fac <- as.double(height)/max(query$layers$height)
  } else {
    stop("Provide either a ref or height object!")
  }

  query$layers$height <- query$layers$height * fac
  query$hs <- tail(query$layers$height, n = 1)

  if ("changes" %in% names(query)) {
    old_changes <- paste0(query$changes, " -> ")
  } else  {
    old_changes <- ""
    query$changes <- paste0(old_changes, "scaleSnowHeight by ", fac)
  }

  return(list(queryScaled = query, trueHeightFactor = (1/fac)))
}


## --- Flip layers top down ----
#' Flip snow profile layers top down
#' @param x snowprofile or snowprofileLayers object with layers to be flipped
#' @return same object with layers dataframe flipped upside down
#' @note only do that with a specific reason (better, don"t do it!), as all functions with snowprofile objects are
#' designed to have the layers increase in height.
#' @export
flipLayers <- function(x) {

  if (is.snowprofileLayers(x)) lay <- x
  else if (is.snowprofile(x)) lay <- x$layers
  else stop("Object must be a snowprofile or snowprofileLayers object")

  lay <- lay[rev(rownames(lay)), ]

  if (is.snowprofileLayers(x)) return(lay)
  else if (is.snowprofile(x)) {
    x$layers <- lay
    return(x)
  }
}




