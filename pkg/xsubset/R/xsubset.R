##
## File:  xsubset.R
##


## interface for fitted lm regression
##
## Args:
##   object       - (lm)
##   model.return - (logical)
##   x.return     - (logical)
##   y.return     - (logical)
##   ...          - arguments forwarded to 'xsubset.default'
##
## Rval:
##   (see 'xsubset.default')
##
## Extracts model matrix, response, weights and offset, and forwards
## the call to 'xsubset.default'.
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
##   object       - (formula) 'lm''s 'formula' argument
##   data         - (data.frame)
##   row.subset   - (numeric[]) 'lm''s 'subset' argument
##   weights      - (numeric[])
##   na.action    - (function)
##   model.return - (logical)
##   x.return     - (logical)
##   y.return     - (logical)
##   contrasts    - (numeric[])
##   offset       - (numeric[])
##   ...          - arguments forwarded to 'xsubset.lm'
##
## Rval: (list)
##   (see 'xsubset.default')
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
##   object    - (matrix)  model matrix
##   y         - (numeric[])  response variable
##   weights   - (numeric[])
##   offset    - (numeric[])
##   include   - (integer[])  regressors to force in
##   exclude   - (integer[])  regressors to force out
##   size      - (integer[])  subset sizes
##   penalty   - (numeric|character)  AIC penalty
##   tolerance - (numeric[])  tolerance
##   pradius   - (integer)  preordering radius
##   nbest     - (integer)  number of best subsets
##   ...       - ignored
##
## Rval: (list)
##   value - (numeric[]|numeric[,])  deviance
##   which - (logical[,]|logical[,,])  subsets
##   ...
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
  intercept <- x.names[1] == "(Intercept)"
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
    rval$tolerance <- tolerance
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
    ## tolerance
    rval$tolerance <- tolerance[1]
    ## rss
    rval$rss <- array(C_rval$rss, dim = nbest, dimnames = list(1:nbest))
    ## aic
    rval$aic <- array(C_rval$aic, dim = nbest, dimnames = list(1:nbest))
    ## sel
    rval$which <- array(NA, dim = c(nvar, nbest),
                        dimnames = list(x.names, 1:nbest))
    rval$which[which0, ] <- C_rval$which[1:(nvar0 * nbest)]
    ## size
    rval$size <- apply(rval$which, 2, sum)
    ## class
    class(rval) <- c("xsubset2", class(rval))
  }

  # done
  rval
}


## print method for 'xsubset' objects
##
## Args:
##   x   - (xsubset)
##   ... - ignored
##
## Rval: (xsubset) invisible
##
print.xsubset <- function (x, ...)
{
  catln <- function (...) base::cat(..., "\n", sep = "")
  
  catln()
  catln("Call:")
  catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))
  
  val <- as.matrix(c(format(x$nvar),
                     if (x$intercept) "yes" else "no",
                     paste(x$include, collapse = " "),
                     paste(x$exclude, collapse = " "),
                     paste(x$size, collapse = " "),
                     "RSS",
                     paste(x$tolerance, collapse = " "),
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
##   x   - (xsubset2)
##   ... - ignored
##
## Rval: (xsubset2) invisible
##
print.xsubset2 <- function (x, ...)
{
  catln <- function (...) base::cat(..., "\n", sep = "")
  
  catln()
  catln("Call:")
  catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))
  
  val <- as.matrix(c(format(x$nvar),
                     if (x$intercept) "yes" else "no",
                     paste(x$include, collapse = " "),
                     paste(x$exclude, collapse = " "),
                     paste("AIC (k = ", x$penalty, ")", sep = ""),
                     format(x$tolerance),
                     format(x$nbest),
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
## Rval: (xsubset) invisible
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
## Rval: (xsubset2) invisible
##
## Graphical arguments are passed to 'plot.default'.
##
plot.xsubset2 <- function (x, type = "b", main = "Deviance", xlab = NULL,
                           ylab = "", col = "blue", lty = 1, legend = TRUE,
                           ...)
{
  ## aic
  aic <- AIC(x, rank = 1:x$nbest)

  ## sub title
  sub <- paste("AIC (k = ", x$penalty, ")", sep = "")

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
  text(rank, aic, labels = paste("[", x$size, "]", sep = ""),
       adj = c(0, 1.5), cex = 0.8)

  ## legend
  if (legend) {
    legend("topleft", "AIC [Size]", lty = lty, col = col, bty = "n")
  }

  ## done
  invisible(x)
}


## extract variable names
##
## Args:
##   object - (xsubset)
##   size   - (integer[])
##   rank   - (integer)
##   ...    - ignored
##
## Rval: (character[]) variable names
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
## Rval: (character[]) variable names
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
##   ...  - forwarded
##
## Rval: (formula)
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
## Rval: (formula)
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
## Rval: (formula)
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
## Rval: (data.frame)
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
## Rval: (data.frame)
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
## Rval: (data.frame)
##
model.frame.xsubset2 <- function (formula, ...) {
  NextMethod(.Generic, formula)
}


## extract model matrix (base method)
##
## Args:
##   object - (xsubset0)
##   ...    - forwarded
##
## Rval: (matrix)
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
## Rval: (matrix)
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
## Rval: (matrix)
##
model.matrix.xsubset2 <- function (object, rank = 1, ...) {
  NextMethod(.Generic, object, rank = rank)
}


## refit (base method)
##
## Args:
##   object    - (xsubset)
##   ...       - forwarded
##   mask.call - (logical) used internally
##
## Rval: (lm)
##
## Note:  'mask.call'
## Reconstruct 'call' object for pretty printing.
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
## Rval: (lm)
##
## Note:  'mask.call'
## Reconstruct 'call' object for pretty printing.
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
## Rval: (lm)
##
## Note:  'mask.call'
## Reconstruct 'call' object for pretty printing.
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
##   ...    - forwarded
##
## Rval: (numeric[])
##
coef.xsubset <- function (object, size, rank = 1, ...) {
  coef(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract ceofficients
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - forwarded
##
## Rval: (numeric[])
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
##   ...    - forwarded
##
## Rval: (matrix)
##
vcov.xsubset <- function (object, size, rank = 1, ...) {
  vcov(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract variance-covariance matrix
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - forwarded
##
## Rval: (matrix)
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
##   ...    - forwarded
##
## Rval: (numeric[])
##
fitted.xsubset <- function (object, size, rank = 1, ...) {
  fitted(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract fitted values
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - forwarded
##
## Rval: (numeric[])
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
##   ...    - forwarded
##
## Rval: (numeric[])
##
residuals.xsubset <- function (object, size , rank = 1, ...) {
  residuals(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract residuals
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - forwarded
##
## Rval: (numeric[])
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
## Rval: (numeric[])
##
## Returns the deviance for the specified subset size and rank.
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
## Returns the deviance for the specified subset size and rank.
##
deviance.xsubset2 <- function (object, rank = 1, ...) {
  if (is.null(rank)) {
    rank <- 1:object$nbest
  }
  object$rss[rank]
}


## extract log-likelihood
##
## Args:
##   object - (xsubset)
##   ...    - forwarded
##
## Rval: (logLik)
##
## Note:  'size'
## This method returns a vector if 'length(size)' is
## greater than one.
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
##   size   - (integer)
##   rank   - (integer)
##   ...    - ignored
##
## Rval: (logLik)
##
logLik.xsubset <- function (object, size, rank = 1, ...) {
  NextMethod(.Generic, object, size = size, rank = rank, df = size + 1)
}


## extract log-likelihood
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - ignored
##
## Rval: (logLik)
##
logLik.xsubset2 <- function (object, rank = 1, ...) {
  ## done
  NextMethod(.Generic, object, rank = rank, df = object$size[rank] + 1)
}


## compute AIC (base method)
##
## Args:
##   object - (xsubset)
##   ...    - forwarded
##   k      - (integer)
##
## Rval: (numeric|data.frame)
##
## Note:  'size'
## The default behavior is to extract the AIC for
## all selected subset sizes.
##
AIC.xsubset0 <- function (object, ..., k = 2) {
  ll <- logLik(object, ...)
  aic <- AIC(ll, k = k)

  if (length(aic) > 1) {
    aic <- data.frame(df = attr(ll, "df"), AIC = aic)
  }

  aic
}


## compute AIC
##
## Args:
##   object - (xsubset)
##   size   - (integer[])
##   rank   - (integer)
##   ...    - ignored
##   k      - (integer)
##
## Rval: (numeric|data.frame)
##
## Note:  'size'
## The default behavior is to extract the AIC for
## all selected subset sizes.
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
##   k      - (integer)
##
## Rval: (numeric|data.frame)
##
AIC.xsubset2 <- function (object, rank = 1, ..., k = NULL) {
  ## return existing values
  if (is.null(k)) {
    return (object$aic[rank])
  }

  ## done
  NextMethod(.Generic, object, rank = rank, k = k)
}


##
## SUMMARY
##


## summary for 'xsubset' objects
##
## Args:
##   object  - (xsubset)
##   rank    - (integer)
##   penalty - (numeric|character)
##   ...     - ignored
##
## Rval: (summary.xsubset)
##
summary.xsubset <- function (object, rank = 1, penalty = 2, ...) {
  paste <- function (...) base::paste(..., sep = "")

  ## rss
  rss <- deviance(object, rank = rank)

  ## aic
  aic <- AIC(object, size = object$size, rank = rank, k = penalty)
  if (length(object$size) > 1) {
    aic <- aic$AIC
    best <- which.min(aic)
  } else {
    best <- 1
  }

  ## names
  names(aic) <- object$size

  ## update object
  object$sum.rank <- rank
  object$sum.penalty <- penalty
  object$sum.rss <- rss
  object$sum.aic <- aic
  object$sum.best <- best
  class(object) <- c("summary.xsubset", class(object))

  ## done
  object
}



## summary for 'xsubset2' objects
##
## Args:
##   object      - (xsubset2)
##   ...         - ignored
##
## Rval: (summary.xsubset2)
##
summary.xsubset2 <- function (object, ...) {
  paste <- function (...) base::paste(..., sep = "")

  ## rss
  rss <- deviance(object, rank = NULL)

  ## update object
  object$rss <- rss
  class(object) <- c("summary.xsubset2", class(object))

  ## done
  object
}


## print 'xsubset' summary
##
## Args:
##   x      - (summary.xsubset)
##   digits - (integer)
##   ...    - ignored
##
## Rval: (summary.xsubset) invisible
##
print.summary.xsubset <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
  catln <- function (..., sep = "") base::cat(..., "\n", sep = sep)
  paste <- function (..., sep = "") base::paste(..., sep = sep)

  ## call
  catln()
  catln("Call:")
  catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))

  ## rank
  catln()
  catln("Rank: ", paste(x$sum.rank, "."), "/", x$nbest)

  ## variable table
  catln()
  catln("Selected variables (by size):")
  which.x <- ifelse(x$which[, x$sum.rank, ][, x$size, drop = FALSE], "x", "")
  rownames(which.x)[x$include] <- paste("+", rownames(which.x)[x$include])
  rownames(which.x)[x$exclude] <- paste("-", rownames(which.x)[x$exclude])
  colnames(which.x)[x$sum.best] <- paste(colnames(which.x)[x$sum.best], "*")
  print(which.x, quote = FALSE)

  ## fit
  catln()
  catln("Model fit:")
  fit <- rbind(x$sum.aic, x$sum.rss)
  rownames(fit) <- c("AIC", "RSS")
  colnames(fit)[x$sum.best] <- paste(colnames(fit)[x$sum.best], "*")
  print(fit, digits = digits)
  catln()
  catln("AIC: penalty = ", format(x$sum.penalty, digits = digits))
  catln()
  catln("* best subset")
    
  invisible(x)
}


## print 'xsubset2' summary
##
## Args:
##   x      - (summary.xsubset2)
##   digits - (integer)
##   ...    - ignored
##
## Rval: (summary.xsubset2) invisible
##
print.summary.xsubset2 <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
  catln <- function (..., sep = "") base::cat(..., "\n", sep = sep)
  paste <- function (..., sep = "") base::paste(..., sep = sep)

  ## call
  catln()
  catln("Call:")
  catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))

  ## rank
  rank <- 1:x$nbest
  catln()
  catln("Rank: ", paste(paste(rank, "."), collapse = " "), "/", x$nbest)

  ## variables
  catln()
  catln("Selected variables (by rank):")
  which.x <- ifelse(x$which, "x", "")
  rownames(which.x)[x$include] <- paste("+", rownames(which.x)[x$include])
  rownames(which.x)[x$exclude] <- paste("-", rownames(which.x)[x$exclude])
  colnames(which.x) <- paste(colnames(which.x), ".")
  print(which.x, quote = FALSE)

  ## fit
  catln()
  catln("Model fit:")
  fit <- rbind(format(x$aic), format(x$rss), format(x$size[rank]))
  rownames(fit) <- c("AIC", "RSS", "Size")
  colnames(fit) <- paste(colnames(fit), ".")
  print(fit, digits = digits, quote = FALSE)
  catln()
  catln("AIC: penalty = ", x$penalty)
    
  ## done
  invisible(x)
}


## plot 'xsubset' summary
##
## Args:
##   x      - (xsubset)
##   type   - (character)
##   main   - (character)
##   xlab   - (character)
##   ylab   - (character)
##   col    - (integer[]|character[])
##   lty    - (integer)
##   legend - (logical)
##   ...    - forwarded
##
## Rval: (summary.xsubset) invisible
##
## All arguments are passed to 'plot.default'.
##
plot.summary.xsubset <- function (x, type = "b", main = NULL, xlab = NULL,
                                  ylab = "", col = c("blue", "red"), lty = 1,
                                  legend = TRUE, ...)
{
  ## x label
  if (is.null(xlab)) {
    xlab <- "Number of regressors in model"
  }

  ## type
  type <- rep(type, length.out = 2)

  ## main title
  if (is.null(main)) {
    main <- "AIC and residual sum of squares"
  }

  ## sub title
  sub <- paste("AIC (k = ", x$sum.penalty, ")", sep = "")

  ## color
  col <- rep(col, length.out = 2)

  ## line type
  lty <- rep(lty, length.out = 2)

  ## plot rss
  plot(x$size, x$sum.rss, type = type[1], main = main, sub = sub,
       xlab = xlab, ylab = ylab, col = col[1], lty = lty[1], ...)

  ## plot aic
  new.old <- getOption("new")
  par(new = TRUE)
  plot(x$size, x$sum.aic, type = type[2], xlab = "", ylab = "",
       col = col[2], lty = lty[2], axes = FALSE, ...)
  par(new = new.old)

  ## legend
  if (legend) {
    legend("top", c("RSS", "AIC"), lty = lty, col = col, bty = "n")
  }

  ## axes
  axis(4)

  ## done
  invisible(x)
}


## plot 'xsubset2' summary
##
## Args:
##   x      - (xsubset2)
##   type   - (character) plot type
##   main   - (character) main title
##   xlab   - (character) x label
##   ylab   - (character) y label
##   col    - (integer[]|character[]) color
##   lty    - (integer) line type
##   legend - (logical) legend?
##   ...    - forwarded
##
## Rval: (summary.xsubset2) invisible
##
## All arguments are passed to 'plot.default'.
##
plot.summary.xsubset2 <- function (x, type = "b", main = NULL, xlab = NULL,
                                  ylab = "", col = c("blue", "red"), lty = 1,
                                  legend = TRUE, ...)
{
  ## type
  type <- rep(type, length.out = 2)

  ## main title
  if (is.null(main)) {
    main <- paste("AIC and residual sum of squares")
  }

  ## sub title
  sub <- paste("AIC (penalty = ", x$penalty, ")", sep = "")

  ## x label
  if (is.null(xlab)) {
    xlab <- "Rank"
  }

  ## color
  col <- rep(col, length.out = 2)

  ## line type
  lty <- rep(lty, length.out = 2)

  ## rank
  rank <- 1:x$nbest

  ## plot rss
  plot(rank, x$rss, type = type[1], main = main, sub = sub,
       xlab = xlab, ylab = ylab, col = col[1], lty = lty[1], ...)

  ## plot aic
  new.old <- getOption("new")
  par(new = TRUE)
  plot(rank, x$aic, type = type[2], xlab = "", ylab = "",
       col = col[2], lty = lty[2], axes = FALSE, ...)
  par(new = new.old)

  ## legend
  if (legend) {
    legend("top", c("RSS", "AIC"), lty = lty, col = col, bty = "n")
  }

  ## axes
  axis(4)

  ## done
  invisible(x)
}


## summary for 'summary.xsubset' objects
##
## Args:
##   object  - (summary.xsubset)
##   ...     - forwarded
##
## Rval: (summary.xsubset)
##
## Note:
## Utility method to tidy class attribute
##
summary.summary.xsubset <- function (object, ...) {
  ## remove 'summary.xsubste' from classes
  class(object) <- class(object)[-1]

  ## done
  NextMethod(.Generic, object, ...)
}


## summary for 'summary.xsubset2' objects
##
## Args:
##   object  - (summary.xsubset2)
##   ...     - forwarded
##
## Rval: (summary.xsubset2)
##
## Note:
## Utility method to tidy class attribute
##
summary.summary.xsubset2 <- function (object, ...) {
  ## remove 'summary.xsubset' from classes
  class(object) <- class(object)[-1]

  ## done
  NextMethod(.Generic, object, ...)
}


## refit
##
## Args:
##   object    - (xsubset)
##   ...       - ignored
##   mask.call - (logical) used internally
##
## Rval: (lm)
##
## Note:  'mask.call'
## Reconstruct 'call' object for pretty printing.
##
refit.summary.xsubset <- function (object, ..., mask.call = TRUE) {
  best <- object$size[object$sum.best]
  NextMethod("summary.xsubset", object, size = best,
             rank = object$sum.rank, mask.call = mask.call)
}


## refit
##
## Args:
##   object    - (xsubset)
##   ...       - ignored
##   mask.call - (logical) used internally
##
## Rval: (lm)
##
## Note:  'mask.call'
## Reconstruct 'call' object for pretty printing.
##
refit.summary.xsubset2 <- function (object, ..., mask.call = TRUE) {
  NextMethod("summary.xsubset2", object, rank = 1, mask.call = mask.call)
}
