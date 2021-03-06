#' Warp one snow profile onto another one
#'
#' After the DTW alignment of two profiles, the maps between the two profiles can be used to
#' warp one profile onto the other profile. In other words, the layer thicknesses of the warped profile
#' are adjusted to optimally align with the corresponding layers of the other profile.
#'
#' After this procedure, the thickness of some layers can be zero, which leads to the layers disappearing.
#'
#' This function is automatically called in `dtwSP(..., keep.internals = TRUE)` to warp the query profile
#' onto the reference profile.
#'
#' *Whom* to warp: There exist 8 different options, 4 for warping the query onto the ref and 4 for vice versa.
#' The 4 options for warping the query onto the ref are:
#'
#'   - global alignment / partial alignment where entire query is matched to subsequence of ref ("jmin")
#'   - partial alignment where entire ref is matched to subsequence of query ("imin")
#'   - partial top down alignment where entire query is matched to subsequence of ref ("jminTopDown")
#'   - partial top down alignment where entire ref is matched to subsequence of query ("iminTopDown")
#'
#' For the other case, warping the ref onto the query, only the equivalent of the first option is implemented.

#'
#' @param alignment DTW alignment object from [dtwSP] containing the two profiles (i.e., called `dtwSP(..., keep.internals = TRUE)`)
#' @param whom whom to warp? "query" (= "jmin"), "imin", "queryTopDown" (= "jminTopDown"), "iminTopDown", "ref";
#' if 'NA' the routine determines that itself from the structure of the alignment object. (see Details)
#'
#' @return Returns the input alignment object including the element alignment$queryWarped (or $referenceWarped),
#' which are the warped snow profiles. The class of the alignment object is altered to "dtwSP", but still inherits "dtw".
#'
#' @author fherla
#'
#' @examples
#'
#' ## first align profiles
#' alignment <- dtwSP(SPpairs$A_modeled, SPpairs$A_manual, open.end = FALSE)
#'
#' ## warp reference profile onto query profile:
#' refWarped <- warpSP(alignment, whom = "ref")$referenceWarped
#' opar <- par(no.readonly =TRUE)
#' par(mfrow = c(1, 2))
#' plot(alignment$query, main = "query")
#' plot(refWarped, main = "warped reference")
#' par(opar)
#'
#' @export
warpSP <- function(alignment, whom = NA) {

  ## find out whom to warp from the alignment:
  if (is.na(whom)) {
    whom <- alignment$openEndType
  }
  if (alignment$direction == "topDown") whom <- paste0(whom, "TopDown")

  q <- alignment$query
  r <- alignment$reference

  ## --- Warp QUERY ----------------------------------------------------------------------------------------------------
  if (whom %in% c("query", "jmin")) {
    ## warping approach analogously to example(dtw) using dtw-object index vectors:
    qmod <- q[names(q) != "layers"]
    qmodHeight <- r$layers$height[alignment$index2]
    qmod$layers <- snowprofileLayers(gtype = q$layers$gtype[alignment$index1],
                                     hardness = q$layers$hardness[alignment$index1],
                                     height = qmodHeight,
                                     formatTarget = "thickness")
    class(qmod) <- class(q)
    qmod <- rmZeroThicknessLayers(qmod)
    alignment$queryWarped <- qmod

  } else if (whom == "imin") {
    ## OE warp that matched the subsequence of the query profile:
    ## i.e. warp query, but stack non-matched layers on top
    qmod <- q[names(q) != "layers"]
    qmodHeight <- q$layers$height[alignment$index2]
    qmod$layers <- snowprofileLayers(gtype = as.factor(c(as.character(q$layers$gtype[alignment$index1]),
                                                    as.character(q$layers$gtype[-(1:alignment$imin-1)]))),
                                     hardness = c(q$layers$hardness[alignment$index1],
                                                  q$layers$hardness[-(1:alignment$imin-1) ]),
                                     height = c(qmodHeight,
                                                qmodHeight[length(qmodHeight)] - q$layers$height[alignment$imin] + q$layers$height[-(1:alignment$imin-1)]),
                                     formatTarget = "thickness")
    class(qmod) <- class(q)
    qmod <- rmZeroThicknessLayers(qmod)
    alignment$queryWarped <- qmod

    ## :::::::::::::::::::::::::::::::::::::::::::::::::::
    ## Top Down
    ## :::::::::::::::::::::::::::::::::::::::::::::::::::
  } else if (whom %in% c("queryTopDown", "jminTopDown")) {
    ## need to reverse vectors in top down:
    qmod <- q[names(q) != "layers"]
    qmodHeight <- rev(r$layers$height[alignment$index2])
    qmod$layers <- snowprofileLayers(gtype = rev(q$layers$gtype[alignment$index1]),
                                     hardness = rev(q$layers$hardness[alignment$index1]),
                                     height = qmodHeight,
                                     formatTarget = "thickness")
    class(qmod) <- class(q)
    qmod <- rmZeroThicknessLayers(qmod)

    alignment$queryWarped <- qmod


  } else if (whom == "iminTopDown") {
    ## need to reverse vectors and undo stacking:
    qmod <- q[names(q) != "layers"]
    qmodHeight <- rev(q$layers$height[alignment$index2])
    qmod$layers <- snowprofileLayers(gtype = as.factor(rev(as.character(q$layers$gtype[alignment$index1]))),
                                     hardness = rev(q$layers$hardness[alignment$index1]),
                                     height = qmodHeight,
                                     formatTarget = "thickness")
    class(qmod) <- class(q)
    qmod <- rmZeroThicknessLayers(qmod)

    alignment$queryWarped <- qmod

    ## --- Warp REF ----------------------------------------------------------------------------------------------------
  } else if (whom %in% c("ref", "reference")) {
    warning("warpSP: warping the reference object is currently in beta stage, check whether output makes sense!")
    ## analogous to "query" (but with indices swapped):
    rmod <- r[names(r) != "layers"]
    rmodHeight <- q$layers$height[alignment$index1]
    rmod$layers <- snowprofileLayers(gtype = r$layers$gtype[alignment$index2],
                                     hardness = r$layers$hardness[alignment$index2],
                                     height = rmodHeight,
                                     formatTarget = "thickness")
    class(rmod) <- class(r)
    rmod <- rmZeroThicknessLayers(rmod)

    alignment$referenceWarped <- rmod


  } else {
    stop(paste0("warpSP: Don't know whom to warp! -> ", whom, "?"))
  }

  ## --- RETURN ----
  if (!inherits(alignment, "dtwSP")) class(alignment) <- append("dtwSP", class(alignment))

  ## clean up data.frame rownames:
  if ("queryWarped" %in% names(alignment)) rownames(alignment$queryWarped$layers) <- seq(nrow(alignment$queryWarped$layers))
  else if ("referenceWarped" %in% names(alignment)) rownames(alignment$referenceWarped$layers) <- seq(nrow(alignment$referenceWarped$layers))

  return(alignment)
}
