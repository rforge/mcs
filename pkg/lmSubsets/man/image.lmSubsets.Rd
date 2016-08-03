\name{image.lmSubsets}
\alias{image.lmSubsets}

\title{Heatmaps for Visualizing Subset Selections}

\description{
  Visualizations of best subset selections carried out
  via \code{\link{lmSubsets}}.
}

\usage{
\method{image}{lmSubsets}(x, nbest = 1, size = NULL, which = NULL,
  xlab = "", ylab = NULL, main = "Subset selection",
  xaxs = "i", yaxs = "i", cex = 0.9,
  col = gray.colors(2), rcol = 2, rect = "BIC", uline = rect,
  srt = 45, ...)

}

\arguments{
  \item{x}{an object of class \code{\link{lmSubsets}}.}
  \item{nbest}{Number of best subset solution.}
  \item{size}{Subset sizes to be plotted (default: all).}
  \item{which}{Regressors to be plotted (default: all).}
  \item{xlab, ylab, main}{Axis and main titles.}
  \item{xaxs, yaxs}{The style of axis interval calculation to be used (see 
    \code{\link{par}} for details).}
  \item{cex}{Expansion factor for labels.}
  \item{col}{A vector of colors of length two used for coloring the
    selected and non-selected regressors, respectively.}
  \item{rcol}{Color for the selected \code{rect}angle(s) (if any).}
  \item{rect, uline}{Logical, character, or numeric for specifying
    which subset size(s) should be highlighted by a rectangle
    and underlining, respectively.}
  \item{srt}{Angle of the x-axis labels.}
  \item{\dots}{Further arguments passed to \code{\link[graphics]{plot}}.}
}

\seealso{\code{\link{lmSubsets}}}

\examples{
## all-subsets regression
data("AirPollution", package = "lmSubsets")
lms <- lmSubsets(mortality ~ ., data = AirPollution)

## visualize the selected subsets for all sizes and
## highlight best BIC solution (default)
image(lms)

## variations
## only some subset sizes and no-highlighting
image(lms, size = 3:10, rect = FALSE)
## highlighting of multiple sizes
image(lms, rect = c(5, 10), rcol = c(2, 4))
## rectangle for subset size 3 but not underlining (and vice versa)
image(lms, rect = 3, uline = FALSE)
image(lms, rect = FALSE, uline = 3)
}

\keyword{regression}
