##
## File:  xsubset.R
##


xsubset <- function (object, ...)
  UseMethod("xsubset")


refit <- function (object, ...)
  UseMethod("refit")



################
## GENERATORS ##
################


## interface for fitted lm regression
##
## Args:
##   object       - (lm)
##   model.return - (logical)
##   x.return     - (logical)
##   y.return     - (logical)
##   ...          - forwarded to 'xsubset.default'
##
## Rval:  (xsubset0)
##   See 'xsubset.default'.
##
xsubset.lm <- function (object, model.return = TRUE, x.return = FALSE,
                        y.return = FALSE, ...)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## terms and model frame
  mt <- terms(object)
  mf <- model.frame(object)

  ## model matrix
  if (is.null(x <- object[["x"]])) {
    x <- model.matrix(mt, mf)
  }

  ## model response
  if (is.null(y <- object[["y"]])) {
    y <- model.response(mf)
  }
  y.name <- deparse(attr(mt, "variables")[[attr(mt, "response") + 1]])

  ## weights and offset
  w <- weights(object)
  o <- object$offset

  ## forward call
  rval <- xsubset(x, y, weights = w, offset = o, ...)

  ## update return value
  rval$call <- call
  rval$formula <- formula(object)
  rval$na.action <- object$na.action
  if (model.return) rval$model <- mf
  if (x.return) rval$x <- x
  if (y.return) rval$y <- y
  rval$contrasts <- object$contrasts
  rval$terms <- mt
  rval$xlevels <- object$xlevels
  rval$y.name <- y.name

  ## done
  rval
}


## standard formula interface
##
## Args:
##   formula      - (formula)
##   data         - (data.frame)
##   row.subset   - (numeric[])
##   weights      - (numeric[])
##   na.action    - (function)
##   model.return - (logical)
##   x.return     - (logical)
##   y.return     - (logical)
##   contrasts    - (numeric[])
##   offset       - (numeric[])
##   ...          - forwarded to 'xsubset.default'
##
## Rval:  (xsubset0)
##   See 'xsubset.default'.
##
xsubset.formula <- function (formula, data, row.subset, weights,
                             na.action = na.omit, model.return = TRUE,
                             x.return = FALSE, y.return = FALSE,
                             contrasts = NULL, offset, ...)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1L]] <- as.name("xsubset")

  ## model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$subset <- mf$row.subset;  mf$row.subset <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## terms
  mt <- attr(mf, "terms")

  ## empty model?
  if (is.empty.model(mt)) {
    ## TODO
    stop ("empty model")
  }

  ## model matrix
  x <- model.matrix(mt, mf, contrasts)

  ## response
  y <- model.response(mf)
  y.name <- deparse(attr(mt, "variables")[[attr(mt, "response") + 1]])

  ## weights
  w <- model.weights(mf)

  ## offset
  o <- model.offset(mf)

  ## forward call
  rval <- xsubset(x, y, weights = w, offset = o, ...)

  ## update return value
  rval$call <- call
  rval$formula <- formula
  rval$na.action <- attr(mf, "na.action")
  if (model.return) rval$model <- mf
  if (x.return) rval$x <- x
  if (y.return) rval$y <- y
  rval$contrasts <- attr(x, "contrasts")
  rval$terms <- mt
  rval$xlevels <- .getXlevels(mt, mf)
  rval$y.name <- y.name

  ## done
  rval
}


## default method
##
## Args:
##   object    - (matrix) model matrix
##   y         - (numeric[]) response variable
##   weights   - (numeric[])
##   offset    - (numeric[])
##   include   - (integer[]) regressors to force in
##   exclude   - (integer[]) regressors to force out
##   size      - (integer[]) subset sizes
##   penalty   - (numeric) AIC penalty
##   tolerance - (numeric[]) tolerance
##   pradius   - (integer) preordering radius
##   nbest     - (integer) number of best subsets
##   ...       - ignored
##
## Rval:  (xsubset0)
##   call      - (call)
##   weights   - (numeric[])
##   offset    - (numeric[])
##   nobs      - (integer)
##   nvar      - (integer)
##   x.names   - (character[])
##   y.name    - (character)
##   include   - (integer[])
##   exclude   - (integer[])
##   intercept - (logical)
##   size      - (integer)
##   penalty   - (numeric)
##   tolerance - (numeric[])
##   nbest     - (integer)
##   rss       - (array)
##   aic       - (array)
##   which     - (array)
##   .nodes    - (integer)
##
xsubset.default <- function (object, y, weights = NULL, offset = NULL,
                             include = NULL, exclude = NULL,
                             size = NULL, penalty = 0, tolerance = 0,
                             pradius = NULL, nbest = 1, ...)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## model matrix
  x <- as.matrix(object);  object <- NULL

  ## weights
  w <- if (is.null(weights)) rep(1, NROW(x)) else weights
  wok <- w != 0;  wnz <- w[wok]
  x <- sqrt(wnz) * x[wok, , drop = FALSE]
  y <- sqrt(wnz) * y[wok,   drop = FALSE]

  ## offset
  o <- if (is.null(offset)) rep(0, NROW(x)) else offset[wok]
  y <- y - o
  
  ## dims
  nobs <- NROW(x)
  nvar <- NCOL(x)

  ## variables names
  x.names <- colnames(x)
  if (is.null(x.names)) {
    x.names <- paste("X", 1:nvar, sep = "")
  }

  ## include processing
  if (is.null(include)) {
    include <- numeric(0)
  } else {
    ## character --> numeric
    if (is.character(include)) {
      include <- match(include, x.names)
    }
    ## logical --> numeric
    if (is.logical(include)) {
      include <- which(rep(include, length.out = NCOL(x)))
    }
    ## canonicalize include
    include <- sort(unique(include))
  }
  
  ## exclude processing
  if (is.null(exclude)) {
    exclude <- numeric(0)
  } else {
    ## character --> numeric
    if (is.character(exclude)) {
      exclude <- match(exclude, x.names)
    }
    ## logical --> numeric
    if (is.logical(exclude)) {
      exclude <- which(rep(exclude, length.out = nvar))
    }
    ## canonicalize exclude
    exclude <- sort(unique(exclude))
  }

  ## include/exclude non-overlapping
  if (any(intersect(include, exclude))) {
    warning ("'include' and 'exclude' overlapping; fixing 'exclude'")
    exclude <- setdiff(exclude, include)
  }

  ## handle intercept
  intercept <- isTRUE(all.equal(as.vector(x[, 1]), rep(1, nobs)))
  if (intercept) {
    if (any(include == 1)) {
      ## OK, already selected
    } else if (any(exclude == 1)) {
      ## not selected
      intercept <- FALSE
    } else {
      ## select
      include <- c(1, include)
    }
  }

  ## include/exclude columns
  which0 <- setdiff(1:nvar, exclude)
  x0 <- x[, which0, drop = FALSE]
  nvar0 <- NCOL(x0)

  ## process pradius
  if (is.null(pradius)) {
    pradius <- round(nvar0 / 3)
  }

  ## process penalty
  ## ...

  ## process size
  mark <- length(include)
  if (is.null(size)) {
    size <- (mark + 1):nvar0
  }

  ## process tolerance
  tolerance <- rep(tolerance, length.out = nvar0)
  if (penalty == 0) {
    tolerance[-size] <- .Machine$double.xmax
  }

  ## call C function
  C_args <- list(## in
                 nobs      = as.integer(nobs),
                 nvar      = as.integer(nvar0),
                 xy        = as.double(cbind(x0, y)),
                 mark      = as.integer(mark),
                 penalty   = as.double(penalty),
                 tolerance = as.double(tolerance),
                 pradius   = as.integer(pradius),
                 nbest     = as.integer(nbest),
                 ## out
                 rss       = as.double(rep(0, nvar0 * nbest)),
                 aic       = as.double(rep(0, nvar0 * nbest)),
                 which     = as.logical(rep(0, nvar0 * nbest * nvar0)),
                 nodes     = integer(1))
  C_rval <- do.call(".C", c(name = "R_select", C_args))

  ## return value
  rval <- list(call      = call,
               weights   = weights,
               offset    = offset,
               nobs      = nobs,
               nvar      = nvar,
               x.names   = x.names,
               y.name    = NULL,
               include   = include,
               exclude   = exclude,
               intercept = intercept,
               size      = NULL,
               penalty   = penalty,
               tolerance = NULL,
               nbest     = nbest,
               rss       = NULL,
               aic       = NULL,
               which     = NULL,
               .nodes = C_rval$nodes)
  class(rval) <- "xsubset0"

  ## extract value & subsets
  if (penalty == 0) {
    ## size
    rval$size <- size
    ## tolerance
    rval$tolerance <- array(tolerance, dim = nvar0, dimnames = list(1:nvar0))
    ## rss
    rval$rss <- array(NA, dim = c(nbest, nvar),
                      dimnames = list(1:nbest, 1:nvar))
    rval$rss[, 1:nvar0] <- C_rval$rss
    rval$rss[, -size] <- NA
    ## which
    rval$which <- array(NA, dim = c(nvar, nbest, nvar),
                        dimnames = list(x.names, 1:nbest, 1:nvar))
    rval$which[which0, , 1:nvar0] <- C_rval$which
    rval$which[, , -size] <- NA
    ## pretty tolerance
    rval$tolerance[rval$tolerance == .Machine$double.xmax] <- +Inf
    ## FIXME: enhance C code to flag uninitialized entries
    rval$rss[rval$rss == .Machine$double.xmax] <- NA
    ## class
    class(rval) <- c("xsubset", class(rval))
  } else {
    cnames <- paste(1:nbest, ".", sep = "")
    ## tolerance
    rval$tolerance <- array(rep(tolerance[1], nbest), dim = nbest,
                            dimnames = list(cnames))
    ## rss
    rval$rss <- array(C_rval$rss, dim = nbest, dimnames = list(cnames))
    ## aic
    rval$aic <- array(C_rval$aic, dim = nbest, dimnames = list(cnames))
    ## sel
    rval$which <- array(NA, dim = c(nvar, nbest),
                        dimnames = list(x.names, cnames))
    rval$which[which0, ] <- C_rval$which[1:(nvar0 * nbest)]
    ## size
    rval$size <- array(apply(rval$which, 2, sum), dim = nbest,
                       dimnames = list(cnames))
    ## class
    class(rval) <- c("xsubset2", class(rval))
  }

  # done
  rval
}


## print method for 'xsubset' objects
##
## Args:
##   x      - (xsubset)
##   digits - (integer)
##   ...    - ignored
##
## Rval:  (xsubset) invisible
##
print.xsubset <- function (x, digits = NULL, ...)
{
  catln <- function (...) base::cat(..., "\n", sep = "")

  if (is.null(digits)) {
    digits <- max(3, getOption("digits") - 3)
  }
  
  catln()
  catln("Call:")
  catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))
  
  val <- as.matrix(c(format(x$nvar),
                     if (x$intercept) "yes" else "no",
                     paste(x$include, collapse = " "),
                     paste(x$exclude, collapse = " "),
                     paste(x$size, collapse = " "),
                     "RSS",
                     paste(format(x$tolerance, digits = digits, trim = TRUE),
                           collapse = " "),
                     format(x$nbest)),
                     ncol = 1)
  colnames(val) <- ""
  rownames(val) <- c("Total regressors:", "Intercept:", "Include:", "Exclude:",
                     "Size:", "Criterion:", "Tolerance:", "N best:")

  print(val, quote = FALSE)
  
  invisible(x)
}


## print method for 'xsubset2' objects
##
## Args:
##   x      - (xsubset2)
##   digits - (integer)
##   ...    - ignored
##
## Rval:  (xsubset2) invisible
##
print.xsubset2 <- function (x, digits = NULL, ...)
{
  catln <- function (...) base::cat(..., "\n", sep = "")

  if (is.null(digits)) {
    digits <- max(3, getOption("digits") - 3)
  }
  
  catln()
  catln("Call:")
  catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))
  
  val <- as.matrix(c(format(x$nvar),
                     if (x$intercept) "yes" else "no",
                     paste(x$include, collapse = " "),
                     paste(x$exclude, collapse = " "),
                     paste("AIC (k = ", format(x$penalty, digits = digits),
                           ")", sep = ""),
                     paste(format(x$tolerance, digits = digits),
                           collapse = " "),
                     x$nbest,
                     paste(x$size, collapse = " ")),
                     ncol = 1)
  colnames(val) <- ""
  rownames(val) <- c("Total regressors:", "Intercept:", "Include:", "Exclude:",
                     "Criterion:", "Tolerance:", "N best:", "Size*:")

  print(val, quote = FALSE)

  catln()
  catln("* computed sizes")
  
  invisible(x)
}


## plot method for 'xsubset' objects
##
## Args:
##   x      - (xsubset)
##   rank   - (integer[])
##   type   - (character)
##   main   - (character)
##   xlab   - (character)
##   ylab   - (character)
##   col    - (integer[]|character[])
##   lty    - (integer)
##   legend - (logical)
##   ...    - passed ot 'plot.default'
##
## Rval:  (xsubset) invisible
##
## Graphical arguments are passed to 'plot.default'.
##
plot.xsubset <- function (x, rank = 1, type = "b", main = "Deviance",
                          xlab = NULL, ylab = "", col = "blue", lty = 1,
                          legend = TRUE, ...)
{
  ## rss
  rss <- deviance(x, rank = rank)

  ## sub title
  sub <- paste("RSS (rank: ", rank, "./", x$nbest, ")", sep = "")

  ## xlab
  if (is.null(xlab)) {
    xlab <- "Number of regressors in model"
  }

  ## plot
  plot(x$size, rss, type = type, main = main, sub = sub,
       xlab = xlab, ylab = ylab, col = col, lty = lty, ...)

  ## legend
  if (legend) {
    legend("topright", "RSS", col = col, lty = lty, bty = "n")
  }

  ## done
  invisible(x)
}


## plot method for 'xsubset2' objects
##
## Args:
##   x      - (xsubset2)
##   type   - (character)
##   main   - (character)
##   xlab   - (character)
##   ylab   - (character)
##   col    - (integer[]|character[])
##   lty    - (integer)
##   legend - (logical)
##   ...    - passed ot 'plot.default'
##
## Rval:  (xsubset2) invisible
##
## Graphical arguments are passed to 'plot.default'.
##
plot.xsubset2 <- function (x, type = "b", main = "Deviance", xlab = NULL,
                           ylab = "", col = "blue", lty = 1, legend = TRUE,
                           ...)
{
  digits <- max(3, getOption("digits") - 3)

  ## aic
  aic <- AIC(x, rank = 1:x$nbest)$AIC

  ## sub title
  sub <- paste("AIC (k = ", format(x$penalty, digits = digits),
               ")", sep = "")

  ## x label
  if (is.null(xlab)) {
    xlab <- "Rank"
  }

  ## rank
  rank <- 1:x$nbest

  ## plot
  plot(rank, aic, type = type, main = main, sub = sub,
       xlab = xlab, ylab = ylab, col = col, lty = lty, ...)

  ## labels (subset sizes)
  text(rank, aic, labels = paste("(", x$size, ")", sep = ""),
       adj = c(0, 1.5), cex = 0.8)

  ## legend
  if (legend) {
    legend("topleft", "AIC (size)", lty = lty, col = col, bty = "n")
  }

  ## done
  invisible(x)
}



#############
## METHODS ##
#############


## extract variable names
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - ignored
##
## Rval:  (character[])
##   variable names
##
variable.names.xsubset <- function (object, size, rank = 1, ...) {
  ## compute indices
  wi <- object$which[, rank, size]

  ## done
  object$x.names[wi]
}


## extract variable names
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - ignored
##
## Rval:  (character[])
##   variable names
##
variable.names.xsubset2 <- function (object, rank = 1, ...) {
  ## compute indices
  wi <- object$which[, rank]

  ## done
  object$x.names[wi]
}


## extract formula (base method)
##
## Args:
##   x    - (xsubset0)
##   ...  - forwarded to 'variable.names'
##
## Rval:  (formula)
##
formula.xsubset0 <- function (x, ...) {
  ## formula interface?
  if (is.null(x$formula)) {
    return (NULL)
  }

  ## extract variable names
  y.name <- x$y.name
  x.names <- variable.names(x, ...)

  ## handle intercept
  if (x$intercept) {
    x.names <- x.names[-1]
  }

  ## construct formula
  f <- paste(y.name, "~")
  f <- paste(f, paste(x.names, collapse = "+"))
  f <- paste(f, if (x$intercept) "+1" else "-1")

  ## done
  as.formula(f, environment(x))
}


## extract formula
##
## Args:
##   x    - (xsubset)
##   size - (integer)
##   rank - (integer)
##   ...  - ignored
##
## Rval:  (formula)
##
formula.xsubset <- function (x, size, rank = 1, ...) {
  NextMethod(.Generic, x, size = size, rank = rank)
}


## extract formula
##
## Args:
##   x    - (xsubset2)
##   rank - (integer)
##   ...  - ignored
##
## Rval:  (formula)
##
formula.xsubset2 <- function (x, rank = 1, ...) {
  NextMethod(.Generic, x, rank = rank)
}


## extract model frame (base method)
##
## Args:
##   formula - (xsubset0)
##   ...     - ignored
##
## Rval:  (data.frame)
##
model.frame.xsubset0 <- function (formula, ...) {
  ## formula interface?
  if (is.null(formula$formula)) {
    return (NULL)
  }

  ## did we keep the model frame?
  if (!is.null(mf <- formula$model)) {
    return (mf)
  }

  ## forward call
  NextMethod()
}


## extract model frame
##
## Args:
##   formula - (xsubset)
##   ...     - ignored
##
## Rval:  (data.frame)
##
model.frame.xsubset <- function (formula, ...) {
  NextMethod(.Generic, formula)
}


## extract model frame
##
## Args:
##   formula - (xsubset2)
##   ...     - ignored
##
## Rval:  (data.frame)
##
model.frame.xsubset2 <- function (formula, ...) {
  NextMethod(.Generic, formula)
}


## extract model matrix (base method)
##
## Args:
##   object - (xsubset0)
##   ...    - forwarded to 'variable.names'
##
## Rval:  (matrix)
##
model.matrix.xsubset0 <- function (object, ...) {
  ## formula interface?
  if (is.null(object$formula)) {
    return (NULL)
  }

  ## did we keep the model matrix?
  if (is.null(x <- object[["x"]])) {
    x <- model.matrix(terms(object), model.frame(object),
                      contrasts = object$contrasts)
  }

  ## extract names
  x.names <- variable.names(object, ...)

  ## done
  x[, x.names]
}


## extract model matrix
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - ignored
##
## Rval:  (matrix)
##
model.matrix.xsubset <- function (object, size, rank = 1, ...) {
  NextMethod(.Generic, object, size = size, rank = rank)
}


## extract model matrix
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - ignored
##
## Rval:  (matrix)
##
model.matrix.xsubset2 <- function (object, rank = 1, ...) {
  NextMethod(.Generic, object, rank = rank)
}


## refit (base method)
##
## Args:
##   object    - (xsubset0)
##   ...       - forwarded to 'formula'
##   mask.call - (logical) used internally
##
## Rval:  (lm)
##
## Note:  'mask.call'
##   Reconstruct 'call' object for pretty printing.
##
refit.xsubset0 <- function (object, ..., mask.call = TRUE) {
  ## formula interface?
  if (is.null(object$terms)) {
    return (NULL)
  }

  ## extract formula
  f <- formula(object, ...)

  ## weights and offset
  w <- object$weights
  o <- object$offset

  ## prettify call?
  if (mask.call) {
    ## extract call to object
    cl <- match.call(expand.dots = FALSE)
    obj.call <- cl$object

    ## construct call to model frame
    mf.call <- call("model.frame")
    mf.call$formula <- obj.call

    ## construct call to lm
    lm.call <- call("lm")
    lm.call$formula <- f
    lm.call$data <- mf.call
    lm.call$weights <-w
    lm.call$offset <- o

    ## done
    eval(lm.call, parent.frame())
  } else {
    args <- list(formula = f,
                 data = model.frame(object),
                 weights = w,
                 offset = o)
    do.call("lm", args)
  }
}


## refit
##
## Args:
##   object    - (xsubset)
##   size      - (integer)
##   rank      - (integer)
##   ...       - ignored
##   mask.call - (logical) used internally
##
## Rval:  (lm)
##
refit.xsubset <- function (object, size, rank = 1, ..., mask.call = TRUE) {
  NextMethod(.Generic, object, size = size, rank = rank, mask.call = mask.call)
}


## refit
##
## Args:
##   object    - (xsubset2)
##   rank      - (integer)
##   ...       - ignored
##   mask.call - (logical) used internally
##
## Rval:  (lm)
##
refit.xsubset2 <- function (object, rank = 1, ..., mask.call = TRUE) {
  NextMethod(.Generic, object, rank = rank, mask.call = mask.call)
}


## extract ceofficients
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - forwarded to 'coef.lm'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   This method refits the 'xsubset' object and calls 'coef' on the
##   obtained 'lm' object.  In some circumstances it might be more
##   efficient to call 'refit' directly and execute 'coef' on the
##   obtained 'lm' object.
##
coef.xsubset <- function (object, size, rank = 1, ...) {
  coef(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract ceofficients
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - forwarded to 'coef.lm'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.xsubset'.
##
coef.xsubset2 <- function (object, rank = 1, ...) {
  coef(refit(object, rank = rank, mask.call = FALSE), ...)
}


## extract variance-covariance matrix
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - forwarded to 'vcov.lm'
##
## Rval:  (matrix)
##
## Note: 'refit'
##   See 'coef.xsubset'.
##
vcov.xsubset <- function (object, size, rank = 1, ...) {
  vcov(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract variance-covariance matrix
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - forwarded to 'vcov.lm'
##
## Rval:  (matrix)
##
## Note: 'refit'
##   See 'coef.xsubset'.
##
vcov.xsubset2 <- function (object, rank = 1, ...) {
  vcov(refit(object, rank = rank, mask.call = FALSE), ...)
}


## extract fitted values
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - forwarded to 'fitted.lm'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.xsubset'.
##
fitted.xsubset <- function (object, size, rank = 1, ...) {
  fitted(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract fitted values
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - forwarded to 'fitted.lm'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.xsubset'.
##
fitted.xsubset2 <- function (object, rank = 1, ...) {
  fitted(refit(object, rank = rank, mask.call = FALSE), ...)
}


## extract residuals
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - forwarded to 'residuals.lm'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.xsubset'.
##
residuals.xsubset <- function (object, size , rank = 1, ...) {
  residuals(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract residuals
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - forwarded to 'residuals.lm'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.xsubset'.
##
residuals.xsubset2 <- function (object, rank = 1, ...) {
  residuals(refit(object, rank = rank, .mask.call = FALSE), ...)
}


## extract deviance
##
## Args:
##   object - (xsubset)
##   size   - (integer[])
##   rank   - (integer[])
##   ...    - ignored
##
## Rval:  (numeric[])
##
## Returns the RSS for the specified subset(s).
##
deviance.xsubset <- function (object, size = NULL, rank = 1, ...) {
  if (is.null(size)) {
    size <- object$size
  }
  object$rss[rank, size]
}


## extract deviance
##
## Args:
##   object - (xsubset)
##   rank   - (integer[])
##   ...    - ignored
##
## Rval: (numeric[])
##
## Returns the RSS for the specified subset(s).
##
deviance.xsubset2 <- function (object, rank = 1, ...) {
  if (is.null(rank)) {
    rank <- 1:object$nbest
  }
  object$rss[rank]
}


## extract log-likelihood (base method)
##
## Args:
##   object - (xsubset0)
##   ...    - forwarded to 'deviance'
##   df     - (integer[]) degrees of freedom
##
## Rval: (logLik)
##
## Note:
##   Can handle multiple objects.
##
logLik.xsubset0 <- function (object, ..., df) {
  ## handle weights
  if(is.null(w <- object$weights)) {
    nobs <- object$nobs 
    sw <- 0
  } else {
    wnz <- w[w > 0]
    nobs <- sum(wnz)
    sw <- sum(log(wnz))
  }

  ## extract rss
  rss <- deviance(object, ...)

  ## done
  structure(0.5 * (sw - nobs * (log(2 * pi) + 1 - log(nobs) + log(rss))),
    df = df, class = "logLik")
}


## extract log-likelihood
##
## Args:
##   object - (xsubset)
##   size   - (integer[])
##   rank   - (integer[])
##   ...    - ignored
##
## Rval:  (logLik)
##
## Note:  'size'
##   This method returns a numeric vector of length
##   'length(rank) * length(size)'.
##
logLik.xsubset <- function (object, size = NULL, rank = 1, ...) {
  ## size
  if (is.null(size)) {
    size <- object$size
  }
  ## degrees of freedom
  df <- rep(size + 1, each = length(rank))
  ## done
  NextMethod(.Generic, object, size = size, rank = rank, df = df)
}


## extract log-likelihood
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - ignored
##
## Rval:  (logLik)
##
## Note:  'rank'
##   This method returns a numeric vector of the same length as
##   'rank'.
##
logLik.xsubset2 <- function (object, rank = 1, ...) {
  ## done
  NextMethod(.Generic, object, rank = rank, df = object$size[rank] + 1)
}


## compute AIC (base method)
##
## Args:
##   object - (xsubset0)
##   ...    - forwarded to 'logLik'
##   k      - (integer)
##
## Rval: (numeric|data.frame)
##
AIC.xsubset0 <- function (object, ..., k = 2) {
  ## extract log-likelihoods
  ll <- logLik(object, ...)
  ## compute AICs
  aic <- AIC(ll, k = k)
  ## data frame?
  if (length(aic) > 1) {
    aic <- data.frame(df = attr(ll, "df"), AIC = aic)
  }
  ## done
  aic
}


## compute AIC
##
## Args:
##   object - (xsubset)
##   size   - (integer[])
##   rank   - (integer[])
##   ...    - ignored
##   k      - (integer)
##
## Rval: (numeric|data.frame)
##
## Note:  'size'
##   The default behavior is to extract the AIC for all sizes.
##
AIC.xsubset <- function (object, size = NULL, rank = 1, ..., k = 2) {
  ## default size
  if (is.null(size)) {
    size <- object$size
  }

  #done
  NextMethod(.Generic, object, size = size, rank = rank, k = k)
}


## compute AIC
##
## Args:
##   object - (xsubset2)
##   rank   - (integer[])
##   ...    - ignored
##   k      - (integer) penalty
##
## Rval: (numeric|data.frame)
##
## Note:  'k'
##   If 'k' is null, return existing values; else compute AIC from
##   scratch for the specified penalty and rank(s).
##
AIC.xsubset2 <- function (object, rank = 1, ..., k = NULL) {
  ## return existing values
  if (is.null(k)) {
    aic <- object$aic[rank]

    ## if multiple values, turn into data frame
    if (length(aic) > 1) {
      aic <- data.frame(df = object$size[rank] + 1, AIC = aic)
    }

    ## done
    return (aic)
  }

  ## done
  NextMethod(.Generic, object, rank = rank, k = k)
}
