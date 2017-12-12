\name{image}


\alias{image}
\alias{image.lmSubsets}
\alias{image.lmSelect}


\title{Heatmaps for Visualizing Subset Selections}


\description{
  Visualizations of all/best subset selections carried out via
  \code{\link{lmSubsets}}.
}


\usage{

\method{image}{lmSubsets}(x, main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
      size = NULL, best = 1, which = NULL, col = gray.colors(2),
      lab = "lab", hilite, hilite_penalty,
      hilite_col = heat.colors, hilite_lab = "lab", pad_size = 3,
      pad_best = 1, pad_which = 3, axis_pos = -4, axis_tck = -4,
      axis_lab = -10, ..., axes = TRUE, ann = par("ann"))

\method{image}{lmSelect}(x, main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
      best = NULL, which = NULL, col = gray.colors(2),
      lab = "lab", hilite, hilite_col = heat.colors,
      hilite_lab = "lab", pad_best = 2, pad_which = 2,
      axis_pos = -4, axis_tck = -4, axis_lab = -10, ...,
      axes = TRUE, ann = par("ann"))
}


\arguments{
  \item{x}{An object of class \code{"lmSubsets"} or \code{"lmSelect"}.}
  \item{main, sub, xlab, ylab}{Main, sub- and axis titles.}
  \item{size, best}{Submodels to be plotted.}
  \item{which}{Regressors to be plotted.}
  \item{col, lab}{Color and label style.}
  \item{hilite, hilite_penalty}{Submodels to be highlighted.}
  \item{hilite_col, hilite_lab}{Highlighting style.}
  \item{pad_size, pad_best, pad_which}{Padding.}
  \item{axis_pos, axis_tck, axis_lab}{Position of axes, tick length, and
    position of labels}
  \item{\dots}{Ignored.}
  \item{axes}{Plot axes.}
  \item{ann}{Annotate plot.}
}


\seealso{\code{\link{lmSubsets}}, \code{\link{lmSelect}}.}


\examples{
## all-subsets regression
data("AirPollution", package = "lmSubsets")
lm_all <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 20)

## heatmap
image(lm_all, best = 1:3)

## highlight 5 best (BIC)
image(lm_all, best = 1:3, hilite = 1:5, hilite_penalty = "BIC")

## highlight all (BIC)
image(lm_all, best = 1:3, hilite = NULL, hilite_penalty = "BIC")

## fun with colors
lm_best <- lmSelect(lm_all)
image(lm_best, col = cbind(rainbow(20), gray.colors(20)))

## axis labels
image(lm_best, lab = c("bold(lab)", "lab"), hilite = 1,
      hilite_col = "red", hilite_lab = "underline(lab)")
}


\keyword{regression}
