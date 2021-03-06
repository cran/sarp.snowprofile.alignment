% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergeIdentLayers.R
\name{mergeIdentLayers}
\alias{mergeIdentLayers}
\title{Merge layers with identical properties}
\usage{
mergeIdentLayers(x, properties = c("hardness", "gtype"))
}
\arguments{
\item{x}{a snowprofile or snowprofileLayers object with \emph{height} grid information}

\item{properties}{a character array of layer properties that are considered when searching for identical layers
(e.g., \code{hardness}, \code{gtype}, ...)}
}
\value{
A new \code{snowprofileLayers} object will be returned with the dimensions \code{height}, \code{hardness}, \code{gtype} and any other
properties given in 'properties'. Depth and thickness information will be auto-calculated. For snowprofile objects, the
field 'changes' will be initialized or extended.
}
\description{
Merge adjacent layers that have identical properties, such as grain type, hardness etc..
}
\examples{

## Merge identical layers based on hardness and grain type:
fewerLayers <- mergeIdentLayers(x = SPpairs$A_modeled, properties = c("hardness", "gtype"))
summary(SPpairs$A_modeled)[, c("hs", "nLayers")]
summary(fewerLayers)[, c("hs", "nLayers")]

## compare profile plots before and after merging (i.e., appear identical!)
opar <- par(no.readonly =TRUE)
par(mfrow = c(1, 2))
plot(SPpairs$A_modeled, main = "original", ylab = "Snow height")
plot(fewerLayers, main = "merged layers", ylab = "Snow height")
par(opar)

}
\author{
fherla
}
