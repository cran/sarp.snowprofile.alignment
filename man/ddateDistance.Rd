% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ddateDistance}
\alias{ddateDistance}
\title{Deposition Date Distance}
\usage{
ddateDistance(
  ddate1,
  ddate2,
  normalizeBy = 5,
  clipWindow = FALSE,
  na.dist = 0.5
)
}
\arguments{
\item{ddate1}{1D array of POSIX dates}

\item{ddate2}{same format and length as ddate1}

\item{normalizeBy}{Numeric scalar to be used for normalization, i.e. the number of days, that defines the distance value of 1}

\item{clipWindow}{Should differences larger than 'normalizeBy' number of days be set to distance 'Infinity'? default FALSE.}

\item{na.dist}{replace NA values with that distance}
}
\value{
An array of length(ddate1) containing the distances according to the configurations.
}
\description{
Calculate the distance (i.e. dissimilarity) between two deposition dates
}
\examples{
## create ddate arrays..
ddate <- as.POSIXct("2019/04/20 12:00", tz = "UTC")
ddate1 <- rep(ddate, 5)
ddate2 <- as.POSIXct(c("2019/04/12 08:00", "2019/04/16 10:00", "2019/04/20 12:00",
                       "2019/04/21 16:00", "2019/04/22 20:00"), tz = "UTC")

## .. and calculate distance:
ddateDistance(ddate1, ddate2, normalizeBy = 5)
}
\author{
fherla
}
