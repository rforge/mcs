\name{image.lmSubsets}
\alias{image.lmSubsets}

\title{Heatmaps for Visualizing Subset Selections}

\description{
  Visualizations of best subset selections carried out
  via \code{\link{lmSubsets}}.
}

\usage{
\method{image}{lmSubsets}(x, best = 1, size = NULL, which = NULL,
  main = "Subset selection", xlab = "", ylab = NULL,
  xaxs = "i", yaxs = "i", cex = 0.9,
  col = gray.colors(2), hilite = "BIC", hilite.hue = 0,
  uline = hilite, srt = 45, gaps = "white", ...)

}

\arguments{
  \item{x}{An object of class \code{\link{lmSubsets}}.}
  \item{best}{Rank of best subset solution.}
  \item{size}{Subset sizes to be plotted (default: all).}
  \item{which}{Regressors to be plotted (default: all).}
  \item{main, xlab, ylab}{Main and axis titles.}
  \item{xaxs, yaxs}{The style of axis interval calculation to be used (see 
    \code{\link{par}} for details).}
  \item{cex}{Expansion factor for labels.}
  \item{col}{A vector of colors of length two used for coloring the
    selected and non-selected regressors, respectively.}
  \item{hilite, uline}{Logical, character, or numeric for specifying
    which subset size(s) should be highlighted or underlined,
    respectively.}
  \item{hilite.hue}{Hue for the selected \code{rect}angle(s) (if any).}
  \item{srt}{Angle of the x-axis labels.}
  \item{gaps}{Color of border between tiles.}
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
## only some subset sizes and no highlighting
image(lms, size = 3:10, hilite = FALSE)
## highlight multiple sizes
image(lms, hilite = c(5, 10), hilite.hue = c(0, 0.65))
## highlight subset size 3 but do not underline (and vice versa)
image(lms, hilite = 3, uline = FALSE)
image(lms, hilite = FALSE, uline = 3)
}

\keyword{regression}
