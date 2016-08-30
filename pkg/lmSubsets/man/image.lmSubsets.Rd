\name{image.lmSubsets}
\alias{image.lmSubsets}

\title{Heatmaps for Visualizing Subset Selections}

\description{
  Visualizations of best subset selections carried out
  via \code{\link{lmSubsets}}.
}

\usage{
\method{image}{lmSubsets}(x, best = 1, size = NULL, which = NULL,
  hilite = "BIC", uline = hilite,
  main = "Subset selection", xlab = "", ylab = NULL,
  col = "gray30", tint = 0.8, gaps = "white",
  hilite.col = "red", hilite.tint = tint,
  xaxs = "i", yaxs = "i", cex = 0.9, srt = 45, ...)
}

\arguments{
  \item{x}{An object of class \code{\link{lmSubsets}}.}
  \item{best}{Rank of best subset solution.}
  \item{size}{Subset sizes to be plotted (default: all).}
  \item{which}{Regressors to be plotted (default: all).}
  \item{hilite, uline}{Logical, character, or numeric specifying
    which subset size(s) should be highlighted or underlined,
    respectively.}
  \item{main, xlab, ylab}{Main and axis titles.}
  \item{col, tint}{Color and tint for selected and non-selected
    regressors, respectively.}
  \item{gaps}{Color of border between tiles.}
  \item{hilite.col, hilite.tint}{Color and tint for highlighted subset sizes.}
  \item{xaxs, yaxs}{The style of axis interval calculation to be used (see 
    \code{\link{par}} for details).}
  \item{cex}{Expansion factor for labels.}
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
## only some subset sizes and no highlighting
image(lms, size = 3:10, hilite = FALSE)
## highlight multiple sizes
image(lms, hilite = c(5, 10), hilite.col = c("green", "blue"))
## highlight subset size 3 but do not underline (and vice versa)
image(lms, hilite = 3, uline = FALSE)
image(lms, hilite = FALSE, uline = 3)
}

\keyword{regression}
