#' Restrict the DTW warping window for snow profiles alignment
#'
#' Given a *quadratic* matrix, this function sets all elements of the matrix that are outside the so-called warping
#' window to `NA`. The warping window is a slanted band of constant width around the main diagonal
#' (i.e., *Sakoe-Chiba*-band), and it's size can be controlled with function arguments.
#'
#' **Note** that the function is designed for cost matrices derived from pairs of snow profiles that have equal numbers
#' of layers (cf., [resampleSPpairs]). The function takes many matrix-like inputs, all of which need to be of the size
#' `(profile.size x profile.size)`, i.e., square matrices of the size 'number of layers'.
#'
#' @param iw matrix of integers indicating their row number (cf., `?row`)
#' @param jw matrix of integers indicating their column number (cf., `?col`)
#' @param iheight matrix of query height filled into the columns of the matrix
#' @param jheight matrix of ref height filled into the rows of the matrix
#' @param iddate same as iheight, but containing deposition date information
#' @param jddate same as jheight, but containing deposition date information
#' @param profile.size number of layers in each of the profiles (scalar)
#' @param profile.height snow height of the profiles (scalar)
#' @param window.size percentage of profile.size or profile.height defining the size of the warping window
#' (i.e., the most restrictive of the two will be applied)
#' @param ddate.window.size number of days that exclude layers from the warping window if their deposition dates
#' differ by more than these days
#' @param ... unused---but important to be able to provide other warping functions to [distMatSP]
#'
#' @seealso [dtw::dtwWindowingFunctions]
#'
#' @export
'warpWindowSP' <- function(iw, jw,
                           iheight, jheight,
                           iddate, jddate,
                           profile.size, profile.height,
                           window.size = 0.3, ddate.window.size = Inf, ...) {

  ## initialize boolean matrix to TRUE
  bmat <- matrix(TRUE, profile.size, profile.size)

  ## set elements outside of warp window to FALSE..
  ## based on index (i.e. layer #):
  bmat[abs(iw - jw) > profile.size * window.size] <- FALSE

  ## based on layer height:
  bmat[abs(iheight - jheight) > profile.height * window.size] <- FALSE

  ## based on ddate:
  if (!all(is.na(iddate))) {
    bmat[abs(iddate - jddate) > ddate.window.size] <- FALSE
  }

  return(bmat)
}
