#' Grain Type similarity matrix for DTW alignments
#'
#' Get the relative similarity matrix of grain types as used for snow profile alignments. This similarity matrix
#' considers the formation and metamorphosis of grain types, as well as quirks of the SNOWPACK model. \cr
#' [grainSimilarity_evaluate] is an analogous matrix designed for assessing the similarity between two profiles, which
#' requires considering the resulting avalanche hazard implications of grain types. \cr
#' The domain is `[0, 1]` --- `1` representing identical grain types. The column 'NA' can be used for unknown grain
#' types.
#' @param triag Return a triangular matrix (TRUE, default) or a symmetric matrix (FALSE)
#' @return data.frame, either triangular or symmetric
#' @author fherla
#' @seealso [grainSimilarity_evaluate], [layerWeightingMat]
#' @examples
#'
#' ## "similarity" matrix:
#' simMat <- grainSimilarity_align()
#' print(simMat)
#'
#' ## equivalent "distance" matrix:
#' distMat <- sim2dist(grainSimilarity_align())
#' print(distMat)
#'
#' @export
grainSimilarity_align <- function(triag = TRUE) {
  DF <- data.frame(
    #    c('PP','DF','RG','FC','DH','SH','MF','IF','PPgp', 'PPsd', 'FCxr','MFcr')
    PP = c(1,    0.8, 0.5, 0.2, 0,   0,   0,   0,    1,    1,       0.2,   0),
    DF = c(NA,   1,   0.8, 0.4, 0,   0,   0.1, 0,    0.8,  0.8,     0.4,   0),
    RG = c(NA,   NA,  1,   0.4, 0.1, 0,   0.1, 0,    0.5,  0.5,     0.5,   0),
    FC = c(NA,   NA,  NA,  1,   0.8, 0.6, 0,   0,    0.2,  0.2,     0.8,   0),
    DH = c(NA,   NA,  NA,  NA,  1,   0.9, 0,   0,    0,    0,       0.7,     0),
    SH = c(NA,   NA,  NA,  NA,  NA,  1,   0,   0,    0,    0,       0,     0),
    MF = c(NA,   NA,  NA,  NA,  NA,  NA,  1,   0,    0,    0,       0.1,   0.2),
    IF = c(NA,   NA,  NA,  NA,  NA,  NA,  NA,  1,    0,    0,       0,     0.8),  # new
    PPgp = c(NA, NA,  NA,  NA,  NA,  NA,  NA,  NA,   1,    1,       0.2,   0),  # new (kept same as PP)
    PPsd = c(NA, NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,    1,       0.2,   0),  # new (kept same as PP)
    FCxr = c(NA, NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,    NA,      1,     0),
    MFcr = c(NA, NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,    NA,      NA,    1),  # new
    "na" = c(0.6,0.6, 0.6, 0.5, 0.4, 0.4, 0.5, 0.4, 0.6,   0.6,     0.6,   0.4),
    row.names = c('PP','DF','RG','FC','DH','SH','MF','IF','PPgp','PPsd', 'FCxr','MFcr')
  )
  if (!triag) {
    DF["na",] <- c(DF[, "na"], NA)
    DF[upper.tri(DF)] = t(DF)[upper.tri(DF)]
  }
  return(DF)
}


#' Grain type similarity matrix for evaluation purposes
#'
#' Similar to [grainSimilarity_align], but designed for assessing the similarity between snow profiles based
#' on avalanche hazard relevant characteristics. To be used in combination with [simSP].
#' @param triag Return a triangular matrix (TRUE, default) or a symmetric matrix (FALSE)
#' @return data.frame, either triangular or symmetric
#' @author fherla
#' @examples
#'
#' simMat <- grainSimilarity_evaluate()
#' print(simMat)
#'
#' @export
grainSimilarity_evaluate <- function(triag = TRUE) {
  DF <- data.frame(
    #    c('PP','DF','RG','FC','DH','SH','MF','IF','PPgp', 'PPsd', 'FCxr','MFcr')
    PP = c(1,    0.8, 0.5, 0.2, 0,   0,   0,   0,    1,    1,       0.2,   0),
    DF = c(NA,   1,   0.8, 0.4, 0,   0,   0,   0,    0.8,  0.8,     0.4,   0),
    RG = c(NA,   NA,  1,   0.4, 0.1, 0,   0,   0,    0.5,  0.5,     0.5,   0),
    FC = c(NA,   NA,  NA,  1,   0.4, 0.3, 0,   0,    0.2,  0.2,     0.6,   0),
    DH = c(NA,   NA,  NA,  NA,  1,   0.9,   0,   0,    0,    0,       0.3,     0),
    SH = c(NA,   NA,  NA,  NA,  NA,  1,   0,   0,    0,    0,       0,     0),
    MF = c(NA,   NA,  NA,  NA,  NA,  NA,  1,   0,    0,    0,       0,     0.2),
    IF = c(NA,   NA,  NA,  NA,  NA,  NA,  NA,  1,    0,    0,       0,     0.8),  # new
    PPgp = c(NA, NA,  NA,  NA,  NA,  NA,  NA,  NA,   1,    1,       0.2,   0),  # new (kept same as PP)
    PPsd = c(NA, NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,    1,       0.2,   0),  # new (kept same as PP)
    FCxr = c(NA, NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,    NA,      1,     0),
    MFcr = c(NA, NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,    NA,      NA,    1),  # new
    "na" = c(0.5,0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,   0.5,     0.5,   0.5),
    row.names = c('PP','DF','RG','FC','DH','SH','MF','IF','PPgp','PPsd', 'FCxr','MFcr')
  )
  if (!triag) {
    DF["na",] <- c(DF[, "na"], NA)
    DF[upper.tri(DF)] = t(DF)[upper.tri(DF)]
  }
  return(DF)
}




#' Weighting scheme for preferential layer matching
#'
#' A matrix, of the same form as [grainSimilarity_align], but containing weighting coefficients for preferential layer
#' matching based on the grain types of the layers.
#' @param triag Return a triangular matrix (TRUE, default) or a symmetric matrix (FALSE)
#' @return data.frame, either triangular or symmetric
#' @author fherla
#' @examples
#'
#' weightsMat <- layerWeightingMat()
#' print(weightsMat)
#'
#' @export
layerWeightingMat <- function(triag = TRUE) {
  DF <- data.frame(
    #    c('PP','DF','RG','FC','DH','SH','MF','IF','PPgp', 'PPsd', 'FCxr','MFcr')
    PP = c(1,    1,   1,   1,   1,   1,   1,   1,    1,    1,       1,   1),
    DF = c(NA,   1,   1,   1,   1,   1,   1,   1,    1,    1,       1,   1),
    RG = c(NA,   NA,  1,   1,   1,   1,   1,   1,    1,    1,       1,   1),
    FC = c(NA,   NA,  NA,  1,   0.9, 0.9, 1,   1,    1,    1,       1,   1),
    DH = c(NA,   NA,  NA,  NA,  0,   0,   1,   1,    1,    1,       1,   1),
    SH = c(NA,   NA,  NA,  NA,  NA,  0,   1,   1,    1,    1,       1,   1),
    MF = c(NA,   NA,  NA,  NA,  NA,  NA,  1,   1,    1,    1,       1,   1),
    IF = c(NA,   NA,  NA,  NA,  NA,  NA,  NA,  1,    1,    1,       1,   1),  # new
    PPgp = c(NA, NA,  NA,  NA,  NA,  NA,  NA,  NA,   1,    1,       1,   1),  # new (kept same as PP)
    PPsd = c(NA, NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,    1,       1,   1),  # new (kept same as PP)
    FCxr = c(NA, NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,    NA,      1,   1),
    MFcr = c(NA, NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,    NA,      NA,  0.5),  # new
    "na" = c( 1,  1,   1,   1,   1,   1,   1,   1,   1,     1,      1,   1),
    row.names = c('PP','DF','RG','FC','DH','SH','MF','IF','PPgp','PPsd', 'FCxr','MFcr')
  )
  if (!triag) {
    DF["na",] <- c(DF[, "na"], NA)
    DF[upper.tri(DF)] = t(DF)[upper.tri(DF)]
  }
  return(5*DF)
}
