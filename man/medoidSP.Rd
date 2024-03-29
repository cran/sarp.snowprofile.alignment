% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/medoidSP.R
\name{medoidSP}
\alias{medoidSP}
\title{Find the medoid snow profile among a group of profiles}
\usage{
medoidSP(
  profileList = NULL,
  rescale_resample = TRUE,
  retDistmat = FALSE,
  distmat = NULL,
  verbose = FALSE,
  resamplingRate = 0.5,
  progressbar = requireNamespace("progress", quietly = TRUE),
  ...
)
}
\arguments{
\item{profileList}{List of snowprofile objects}

\item{rescale_resample}{Do you want to uniformly rescale and resample the set of profiles prior to calculating the distance matrix?}

\item{retDistmat}{Do you want to \emph{return} the pairwise distance matrix?}

\item{distmat}{If you have a precalculated distance matrix, provide it here to compute the medoid on it.}

\item{verbose}{print pairwise distance matrix? default FALSE}

\item{resamplingRate}{The resampling rate that is used for the whole set if \code{rescale_resample = TRUE}}

\item{progressbar}{Do you want to print a progress bar with recommended package "progress"?}

\item{...}{arguments passed to \link{distanceSP} and then further to \link{dtwSP}}
}
\value{
If \code{retDistmat = FALSE} return the (named) index of the medoid snow profile, otherwise return a list with the elements
\code{iMedoid} and \code{distmat}.
}
\description{
Find the medoid snowprofile among a group of profiles, based on their pairwise dissimilarity. Either provide a list
of \code{snowprofile} objects, or a precomputed distance matrix. \cr
If you provide a list of profiles the profiles can optionally be rescaled and resampled before the distance matrix
for the medoid calculation is computed. When computing the distance matrix this routine calls \link{distanceSP} for
\emph{every possible pair} of profiles among the group. During that call the profile pair is aligned by \link{dtwSP}
and the aligned pair is evaluated by \link{simSP}.
Note that the number of possible profile pairs grows exponentially with the number of profiles in the group (i.e.,
O(n^2) calls, where n is the number of profiles in the group).
}
\details{
Note that the pairwise distance matrix is modified within the function call to represent a symmetric distance matrix.
That is,, however, not originally the case, since \code{dtwSP(A, B) != dtwSP(B, A)}. The matrix is therefore made symmetric by
setting the similarity between the profiles A and B to \verb{max(\{dtwSP(A, B), dtwSP(B, A)\})}.
}
\examples{
this_example_runs_about_5s <- TRUE
if (!this_example_runs_about_5s) {  # exclude from cran checks

  ## take a list of profiles
  grouplist <- SPgroup[1:4]
  plot(grouplist, SortMethod = 'unsorted', xticklabels = "originalIndices")

  ## calulate medoid profile
  idxMedoid <- medoidSP(grouplist)
  representativeProfile <- grouplist[[idxMedoid]]
  plot(representativeProfile, main = paste0("medoid (i.e., profile ", idxMedoid, ")"))

}
}
\seealso{
\link{reScaleSampleSPx}
}
\author{
fherla
}
