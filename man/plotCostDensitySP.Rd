% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCostDensitySP.R
\name{plotCostDensitySP}
\alias{plotCostDensitySP}
\title{Plot alignment cost density and warping path}
\usage{
plotCostDensitySP(
  alignment,
  localCost = TRUE,
  labelHeight = FALSE,
  marginalPros = TRUE,
  pathCol = "black",
  target = FALSE,
  movingTarget = FALSE,
  tlty = "dotted",
  tlwd = 1.5,
  tcol = "black",
  tcex = 1.5,
  cex.lab = 1,
  xlab = NULL,
  ylab = NULL,
  ...
)
}
\arguments{
\item{alignment}{object from \link{dtwSP}}

\item{localCost}{plot \emph{local} cost matrix, otherwise plot accumulated global cost.}

\item{labelHeight}{plot axes in units of height (cm) or in unitless (i.e., layer index).}

\item{marginalPros}{plot profiles in margins along the axes. default TRUE}

\item{pathCol}{color of warping path}

\item{target}{draw horizontal & vertical lines from matrix cells to corresponding layers in the (marginal) profiles.
Provide either a vector of length 1 (i.e., index of warping path) or length 2 (i.e., x, y coordinates in terms of
layer indices), or a matrix with 2 columns, specifying (x, y) if you desire multiple 'targets'}

\item{movingTarget}{Do you want to draw the warping path only partially, from the origin to the target cross?
Only possible if target cross is given as a scalar! default = FALSE (Useful to create GIF animations of a moving path)}

\item{tlty}{target lty}

\item{tlwd}{target lwd}

\item{tcol}{target col}

\item{tcex}{target cex}

\item{cex.lab}{cex of axis labels (cf. \link{par})}

\item{xlab}{x-axis label to change default labeling}

\item{ylab}{y-axis label to change default labeling}

\item{...}{forwarded to \link{par}}
}
\description{
Plot alignment cost density and warping path, optionally with the two snow profiles plotted in the margins
along the axes.
}
\note{
If you can't see the axis labels, try e.g., \code{par(oma = c(3, 3, 0, 0))} before calling the function. Note, there
seems to be a problem (only sometimes) with the left-hand labels that are for some reason not plotted parallel
to the axis. Also, the routine is not bulletproof with respect to drawing 'targets'. Apologies for any inconveniences!
}
\examples{

## first align profiles:
dtwAlignment <- dtwSP(SPpairs$A_modeled, SPpairs$A_manual, open.end = FALSE)

## then plot cost density:
plotCostDensitySP(dtwAlignment)

## label height instead of layer index, and don't show warping path:
plotCostDensitySP(dtwAlignment, labelHeight = TRUE, pathCol = "transparent")

## draw lines to the cell that corresponds to the DH and SH layers
plotCostDensitySP(dtwAlignment, target = c(191, 208))

## "moving target", i.e., draw warping path only from origin to target:
plotCostDensitySP(dtwAlignment, target = 200, movingTarget = TRUE)
plotCostDensitySP(dtwAlignment, target = 266, movingTarget = TRUE)


## A cool GIF can be created from frames like those
create_GIF <- FALSE
if (create_GIF){
  nPath <- length(dtwAlignment$index1)
  resolution <- 100  # i.e. super low, make value smaller for smoother GIF
  for (k in seq(1, nPath, by = resolution)) {
    plotCostDensitySP(dtwAlignment, target = k, movingTarget = TRUE)
  }
}

}
\author{
fherla
}
