#' Similarity measure between snow profile pairs
#'
#' A simple similarity score --based on avalanche hazard relevant characteristics-- is calculated for two snow profiles
#' that have been aligned onto the same height grid (either through DTW or resampling).
#' If one profile contains more layers than the other one, the layers with a non-matched height
#' represent missing layers and will be treated accordingly. The similarity score is weighted according to grain type
#' classes to accredit the importance of weak layers, new snow and crusts against less important bulk layers.
#' The similarity measure is compatible with top-down alignments and is symmetric with respect to its inputs, i.e.
#' `simSP(P1, P2) == simSP(P2, P1)`.
#' @param ref snowprofile object 1
#' @param qw snowprofile object 2 (matched layers need to be on the same height grid of ref)
#' @param gtype_distMat a distance matrix that stores **distance** information of grain types (*Be careful* to convert
#' similarities, as in [grainSimilarity_evaluate], into dissimilarities with [sim2dist].)
#' @param verbose print similarities of different grain classes to console? default FALSE
#' @param returnDF additionally return the similarities of the grain classes as data.frame (analogously to verbose);
#' the return object then has the fields `$sim` and `$simDF`
#'
#' @return Either a scalar similarity between `[0, 1]` with 1 referring to the two profiles being identical, or
#' (if `returnDF` is TRUE) a list with the elements `$sim` and `$simDF`.
#'
#' @details The grain classes contain the following grain types:
#'   - weak layers (wl): SH and DH
#'   - new snow (pp): PP and DF
#'   - crusts (cr): MFcr and IF
#'   - bulk: the rest (i.e., predominantly RG, FC, FCxr --- MF falls also in here, will maybe be adjusted in future.)
#'
#' @examples
#'
#' ## first align two profiles
#' alignment <- dtwSP(SPpairs$A_modeled, SPpairs$A_manual)
#'
#' ## then assess the similarity of the aligned profiles
#' SIM <- simSP(alignment$queryWarped, alignment$reference, verbose = TRUE)
#'
#' @export
simSP <- function(ref, qw, gtype_distMat = sim2dist(grainSimilarity_evaluate(triag = FALSE)),
                  verbose = FALSE, returnDF = FALSE) {

  ## --- height grid operations ----
  ## get number of layers from both profiles
  refNL <- nrow(ref$layers)
  qwNL <- nrow(qw$layers)
  ## get index vector for layers located at equal heights
  id <- suppressWarnings(which(ref$layers$height == qw$layers$height))  # no warning if they're not of same lengths
  n_id <- length(id)
  ## reverse for top-down alignments:
  if (n_id == 0) {
    ## ref perspective:
    id <- suppressWarnings(which(unique(rev(ref$layers$height)) == unique(rev(qw$layers$height))))
    cref <- length(which(ref$layers$height == ref$layers$height[refNL])) - 1
    id <- rev(refNL + 1 - cref - id)
    ## warped query perspective:
    seq_qw <- suppressWarnings(which(unique(rev(qw$layers$height)) == unique(rev(ref$layers$height))))
    cqw <- length(which(qw$layers$height == qw$layers$height[qwNL])) - 1
    seq_qw <- rev(qwNL + 1 - cqw - seq_qw)
    ## prepend NA rows to get correct index vector:
    if (seq_qw[1] > 1) {
      ref$layers[seq_qw, ] <- ref$layers[id, ]
      ref$layers[seq(seq_qw[1] - 1), ] <- NA
      ref$layers[seq(seq_qw[1] - 1), "height"] <- seq(from = 0, to = qw$layers$height[seq_qw[1]], length.out = seq_qw[1] + 1)[-c(1, seq_qw[1] + 1)]
      ref$layers[seq(seq_qw[1] - 1), "hardness"] <- 0
      id <- seq_qw
    } else if (id[1] > 1) {
      qw$layers[id, ] <- qw$layers[seq_qw, ]
      qw$layers[seq(id[1] - 1), ] <- NA
      qw$layers[seq(id[1] - 1), "height"] <- seq(from = 0, to = ref$layers$height[id[1]], length.out = id[1] + 1)[-c(1, id[1] + 1)]
      qw$layers[seq(id[1] - 1), "hardness"] <- 0
    }
    n_id <- length(id)
  }
  if (n_id > 0) {
    if (any(diff(id) == 0)) {
      warning(paste0("Profiles don't seem to be on the same depth grid!
                     queryWarped: ", paste0(qw$layers$height, collapse = " "), "
                     ref:         ", paste0(ref$layers$height, collapse = " ")))
      verbose <- TRUE  # in case of warning switch to verbose output
    }
  } else {
    warning(paste0("Profiles have no single layer interface at the same height!
                    queryWarped: ", paste0(qw$layers$height, collapse = " "), "
                    ref:         ", paste0(ref$layers$height, collapse = " "), "
                    returning 'NA'.."))
    return(NA)
  }
  ## merge and resample layers that are on same height:
  RES <- resampleSPpairs(qw$layers[id,], ref$layers[id,], mergeBeforeResampling = TRUE, dims = c("gtype", "hardness"))
  rl <- RES$ref
  qwl <- RES$query  # note: RES$query is indeed correct here (naming convention of resampleSPpairs)

  refGrains <- as.character(rl$gtype)
  qwGrains <- as.character(qwl$gtype)
  matchedGrid <- rep(rl$height, times = 2)  # stacked, analogous to matchedDF (further down)
  nGrains <- length(refGrains)

  ## extract non-matched layer information:
  nonMatchedIn <- which(c(refNL, qwNL) > n_id)
  nMI <- length(nonMatchedIn)
  ## merge non-matched layers, but don't issue warning if there is nothing to merge:
  if (nMI == 1) {
    suppressWarnings(missingLayers <- mergeIdentLayers(list(ref, qw)[[nonMatchedIn]]$layers[-(id), ]))
    toDel <- which(missingLayers$height %in% rl$height)
    if (length(toDel) > 0) missingLayers <- missingLayers[-(toDel), ]
    if (nrow(missingLayers) == 0) missingLayers <- NA
  } else if (nMI == 0) {
    missingLayers <- NA
  } else if (nMI == 2) {
    ## distinguish true unmatched layers from zero-thickness layers in case of no depth resampling:
    if (all(ref$layers[-(id), 'height'] %in% rl$height)) nonMatchedIn <- nonMatchedIn[!nonMatchedIn == 1]
    if (all(qw$layers[-(id), 'height'] %in% qwl$height)) nonMatchedIn <- nonMatchedIn[!nonMatchedIn == 2]

    missingLayers <- data.frame()
    for (nM in nonMatchedIn) {
      suppressWarnings(missingLayers <- rbind(missingLayers,
                                              mergeIdentLayers(list(ref, qw)[[nM]]$layers[-(id), ])))
    }
    if (length(nonMatchedIn) > 1) {  # still true unmatched layers in both profiles: that's odd! print to user:
      print("simSP: non-matched layers in both profiles:")
      print(missingLayers)
    }
  }


  ## --- similarity calculations ----
  ## distances and according similarities from various dims:
  dGT <- extractFromScoringMatrix(ScoringFrame = gtype_distMat,
                                  grainType1 = qwGrains,
                                  grainType2 = refGrains)  # vector
  simGT <- sim2dist(dGT)  # vector

  dHHI <- hardnessDistance(qwl$hardness, rl$hardness, normalize = TRUE, absDist = TRUE)  # vector
  simHHI <- sim2dist(dHHI)  # vector

  ## shift similarity to incorporate 'penalty' for bad alignments and 'score' for good ones  (deprecated!)
  ## and combine dimensions
  pen_offset <- 0
  sim <- (simGT * simHHI)[, 1] - pen_offset

  ## separate further evaluation of similarity into grain type categories
  ## (both profiles are important i.e. aim at 'symmetric' similarity score)
  ## (1) unmatched layers get similarity 'indifferent' i.e. no penalty and no score
  if (is.data.frame(missingLayers)) {
    missingDF <- data.frame(grains = as.factor(as.character(missingLayers[, "gtype"])),  # update levels of grain types
                            sim = 0.5)
  } else {
    missingDF <- data.frame()
  }
  ## (2) matched layers:
  ## stack into data.frame, first refGrains, second qwGrains:
  matchedDF <- data.frame(grains = as.factor(c(as.character(refGrains), as.character(qwGrains))),
                          sim = sim)
  if (any(is.na(matchedDF$sim))) warning("simSP: NAs produced in similarity assessment of matched layers. Investigate why!")

  ## combine similarities within the categories
  ## (missing layers get weight factor 2 b/c matched layers come in once for ref and once for qw)
  cat_wls <- c("SH", "DH")
  # cat_fcs <- c("FC")
  cat_cr <- c("MFcr", "IF")
  cat_pps <- c("PP", "DF")
  cat_special <- c(cat_wls, cat_cr, cat_pps)  # will be used exclusively later to get bulk grains
  # cat_special <- c(cat_wls, cat_fcs, cat_cr, cat_pps)  # will be used exclusively later to get bulk grains

  ## eliminate hardness information for WL, MFcr:
  ## matchedDF is stacked with two repetitions of sim (see above)
  iNoHHI <- which(matchedDF$grains %in% c(cat_wls, cat_cr))
  iExc <- which(iNoHHI > nGrains)
  iNoHHI_trans <- iNoHHI
  iNoHHI_trans[iExc] <- iNoHHI[iExc] - nGrains  # translate indices
  matchedDF[iNoHHI, "sim"] <- simGT[iNoHHI_trans, 1] - pen_offset

  ## Calculate specifc weighting scheme for WL, based on how many WLs are apparent in the profiles:
  nWLs <- max(c(length(which(matchedDF$grains[1:nGrains] %in% cat_wls)),
                length(which(matchedDF$grains[(nGrains+1):(2*nGrains)] %in% cat_wls))))
  if (nWLs > 1) {
    ## more than one WL:
    maxMatchedGrid <- max(matchedGrid)
    gridBoundsStep <- maxMatchedGrid/(nWLs+1)
    gridBounds <- seq(gridBoundsStep, maxMatchedGrid+1, length.out = nWLs)
    gridBounds <- matrix(c(0, gridBounds[1:(nWLs-1)], gridBounds), ncol = 2)

    iMatchedWLs <- c(which(matchedDF$grains %in% cat_wls))
    matchedWLs <- matrix(c(matchedGrid[iMatchedWLs], matchedDF[iMatchedWLs, "sim"]), ncol = 2)
    simWLs <- sapply(seq(nWLs), function(i) {
      mean(matchedWLs[which(matchedWLs[, 1] >= gridBounds[i, 1] & matchedWLs[, 1] < gridBounds[i, 2]), 2], na.rm = T)
    })
  } else {
    ## one or none WL:
    simWLs <- mean(matchedDF[matchedDF$grains %in% cat_wls, "sim"])
  }
  ## Calculate specifc weighting scheme for MFcr, based on how many crusts are apparent in the profiles:
  nCRs <- max(c(length(which(matchedDF$grains[1:nGrains] %in% cat_cr)),
                length(which(matchedDF$grains[(nGrains+1):(2*nGrains)] %in% cat_cr))))
  if (nCRs > 1) {
    ## more than one CR:
    maxMatchedGrid <- max(matchedGrid)
    gridBoundsStep <- maxMatchedGrid/(nCRs+1)
    gridBounds <- seq(gridBoundsStep, maxMatchedGrid+1, length.out = nCRs)
    gridBounds <- matrix(c(0, gridBounds[1:(nCRs-1)], gridBounds), ncol = 2)

    iMatchedCRs <- c(which(matchedDF$grains %in% cat_cr))
    matchedCRs <- matrix(c(matchedGrid[iMatchedCRs], matchedDF[iMatchedCRs, "sim"]), ncol = 2)
    simCRs <- sapply(seq(nCRs), function(i) {
      mean(matchedCRs[which(matchedCRs[, 1] >= gridBounds[i, 1] & matchedCRs[, 1] < gridBounds[i, 2]), 2], na.rm = T)
    })
  } else {
    ## one or none CR:
    simCRs <- mean(matchedDF[matchedDF$grains %in% cat_cr, "sim"])
  }



  simDF <- data.frame(
    wl = mean(c(simWLs, suppressWarnings(mean(missingDF[missingDF$grains %in% cat_wls, "sim"]))),
              na.rm = TRUE),
    # fc = mean(c(matchedDF[matchedDF$grains %in% cat_fcs, "sim"],
    #             rep(missingDF[missingDF$grains %in% cat_fcs, "sim"], times = 2))),
    cr = mean(c(simCRs, suppressWarnings(mean(missingDF[missingDF$grains %in% cat_cr, "sim"]))),
              na.rm = TRUE),
    pp = mean(c(matchedDF[matchedDF$grains %in% cat_pps, "sim"],
                   rep(missingDF[missingDF$grains %in% cat_pps, "sim"], times = 1))),
    bulk = mean(c(matchedDF[!matchedDF$grains %in% cat_special, "sim"],
                   rep(missingDF[!missingDF$grains %in% cat_special, "sim"], times = 1)))
  )

  ## shift similarity back to range [0, 1] after computing intra-grainType average:
  simDF <- simDF + pen_offset
  rownames(simDF) <- "sim [0, 1]: "
  ## combine grain type categories to simple similarity:
  simSP <- mean(as.double(simDF[1, ]), na.rm = TRUE)

  ## return NA in case of NULL/problem:
  if (length(simSP) == 0) {
    simSP <- NA
    warning("simSP: problem in calculating similarity score, returning NA")
  }

  simDF <- round((simDF)*100)/100
  if (verbose) {
    print(simDF)
    cat(paste("simple similarity =", round(simSP*1000)/1000, "\n"))
  }

  ifelse(returnDF, return(list(sim = simSP, simDF = simDF)), return(simSP))
}
