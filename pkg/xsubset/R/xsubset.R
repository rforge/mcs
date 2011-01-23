##
## File:  xsubset.R
##


## method for data fame objects
##
## Args:
##   object  - (data.frame) data frame
##   y.which - (integer|character) index of response variable
##   ...     - arguments forwarded to 'xsubset.lm'
##
## Rval: (list)
##   (see 'xsubset.default')
##
## This method simply constructs an 'lm' from the
## object ('lm(object)') and forwards the call to
## 'xsubset.lm'.
##
xsubset.data.frame <- function (object, y.which = 1, ...) {
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## rearrange data frame
  if (is.character(y.which)) {
    y.which <- match(y.which, names(object))
  }
  df <- data.frame(object[y.which], object[-y.which])

  ## forward call
  rval <- xsubset(lm(df), ...)

  ## update return value
  rval$call <- call

  # done
  rval
}


## interface for fitted lm regression
##
## Args:
##   object - (lm)
##   ...    - arguments forwarded to 'xsubset.formula'
##
## Rval:
##   (see 'xsubset.default')
##
## Extracts model frame, weights and offset, and forwards
## the call to 'xsubset.formula'.
##   
xsubset.lm <- function (object, ...)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## forward call
  rval <- xsubset(formula(object), data = model.frame(object),
                  weights = model.weights(object), offset = model.offset(object),
                  ...)

  ## update return value
  rval$call <- call
  
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
##   contrasts    - (numeric[]) ignored
##   offset       - (numeric[])
##   ...          - arguments forwarded to 'xsubset.lm'
##
## Rval: (list)
##   (see 'xsubset.default')
##
xsubset.formula <- function (object, data = NULL, row.subset = NULL,
                             weights = NULL, na.action = na.omit,
                             model.return = TRUE, x.return = FALSE,
                             y.return = FALSE, contrasts = NULL, offset = NULL,
                             ...)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## terms
  t <- terms(object)

  ## data
  if (is.null(data)) {
    data <- environment(t)
  }

  ## model frame
  mf <- model.frame(object, data = data)

  ## model matrix
  x <- model.matrix(t, data = mf)

  ## response
  y <- model.response(mf)

  ## row subset
  s <- row.subset
  if (is.null(s)) {
    s <- rep(TRUE, NROW(x))
  }
  if (is.logical(s)) {
    s <- which(s)
  }
  x <- x[s, , drop = FALSE]
  y <- y[s]

  ## weights
  if (is.null(weights)) {
    w <- rep(1, NROW(x))
  }  else {
    w <- weights[s]
  }

  x0 <- sqrt(w[w != 0]) * x[w != 0, , drop = FALSE]
  y0 <- sqrt(w[w != 0]) * y[w != 0,   drop = FALSE]

  ## offset
  if (is.null(offset)) {
    o <- rep(0, NROW(data))
  } else {
    o <- offset[s]
  }
  y0 <- y0 - o[w != 0]

  ## NAs
  yx0 <- na.action(cbind(y0, x0))

  ## forward call
  rval <- xsubset(yx0[, -1, drop = FALSE],
                  yx0[, 1, drop = FALSE], ...)

  ## update return value
  rval$call <- call
  rval$terms <- t
  if (!is.null(weights)) rval$weights <- w
  if (model.return) rval$model <- mf
  if (x.return) rval$x <- x
  if (y.return) rval$y <- y
  rval$contrasts <- contrasts
  if (!is.null(offset)) rval$offset <- o
  rval$y.name <- deparse(attr(t, "variables")[[attr(t, "response") + 1]])

  ## done
  rval
}


## default method
##
## Args:
##   object    - (matrix)  model matrix
##   y         - (numeric[])  response variable
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
xsubset.default <- function (object, y, include = NULL, exclude = NULL,
                             size = NULL, penalty = 0, tolerance = 0,
                             pradius = NULL, nbest = 1, ...)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## size
  nobs <- NROW(object)
  nvar <- NCOL(object)

  ## variables names
  x.names <- colnames(object)

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
      include <- which(rep(include, length.out = NCOL(object)))
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
  x0 <- object[, which0, drop = FALSE]
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
  tolerance[-size] <- .Machine$double.xmax

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
                 value = as.double(rep(0, nvar0 * nbest)),
                 which = as.logical(rep(0, nvar0 * nbest * nvar0)),
                 nodes = integer(1))
  C_rval <- do.call(".C", c(name = "R_select", C_args))

  ## return value
  rval <- list(call      = call,
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
               value     = NULL,
               which     = NULL,
               .nodes = C_rval$nodes)
  class(rval) <- "xsubset0"

  ## extract value & subsets
  if (penalty == 0) {
    rval$size <- size
    rval$tolerance <- tolerance
    ## rss
    rval$value <- array(NA, dim = c(nbest, nvar),
                        dimnames = list(1:nbest, 1:nvar))
    rval$value[, 1:nvar0] <- C_rval$value
    ## sel
    rval$which <- array(NA, dim = c(nvar, nbest, nvar),
                        dimnames = list(x.names, 1:nbest, 1:nvar))
    rval$which[which0, , 1:nvar0] <- C_rval$which
    ## class
    class(rval) <- c("xsubset", class(rval))
  } else {
    rval$tolerance <- tolerance[1]
    ## aic
    rval$value <- array(C_rval$value, dim = nbest, dimnames = list(1:nbest))
    ## sel
    rval$which <- array(NA, dim = c(nvar, nbest),
                        dimnames = list(x.names, 1:nbest))
    rval$which[which0, ] <- C_rval$which[1:(nvar0 * nbest)]
    ## class
    class(rval) <- c("xsubset2", class(rval))
  }

  ## FIXME: enhance C code to mark uninitialized entries
  rval$value[rval$value == .Machine$double.xmax] <- NA

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
                     format(x$nbest)),
                     ncol = 1)
  colnames(val) <- ""
  rownames(val) <- c("Total regressors:", "Intercept:", "Include:", "Exclude:",
                     "Criterion:", "Tolerance:", "N best:")

  print(val, quote = FALSE)
  
  invisible(x)
}


## plot method for 'xsubset' objects
##
## Args:
##   x      - (xsubset)
##   size   - (integer[])
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
plot.xsubset <- function (x, size = NULL, rank = 1, type = "b",
                          main = "Deviance", xlab = NULL, ylab = "",
                          col = "blue", lty = 1, legend = TRUE, ...)
{
  ## default size
  if (is.null(size)) {
    size <- x$size
  }

  ## rss
  rss <- deviance(x, size = size, rank = rank)

  ## sub title
  sub <- paste("RSS (rank: ", rank, "./", x$nbest, ")", sep = "")

  ## xlab
  if (is.null(xlab)) {
    xlab <- "Number of regressors in model"
  }

  ## plot
  plot(size, rss, type = type, main = main, sub = sub,
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
## Rval: (xsubset2) invisible
##
## Graphical arguments are passed to 'plot.default'.
##
plot.xsubset2 <- function (x, rank = NULL, type = "b", main = "Deviance",
                          xlab = NULL, ylab = "", col = "blue", lty = 1,
                          legend = TRUE, ...)
{
  ## default rank
  if (is.null(rank)) {
    rank <- 1:x$nbest
  }
  
  ## aic
  aic <- deviance(x, rank = rank)

  ## sub title
  sub <- paste("AIC (k = ", x$penalty, ")", sep = "")

  ## x label
  if (is.null(xlab)) {
    xlab <- "Rank"
  }

  ## plot
  plot(rank, aic, type = type, main = main, sub = sub,
       xlab = xlab, ylab = ylab, col = col, lty = lty, ...)

  ## labels (subset sizes)
  size <- paste("[", apply(x$which[, rank], 2, sum), "]", sep = "")
  text(rank, aic, labels = size, adj = c(0, 1.5), cex = 0.8)

  ## legend
  if (legend) {
    legend("topleft", "AIC", lty = lty, col = col, bty = "n")
  }

  ## done
  invisible(x)
}


## extract environment (base method)
##
## Args:
##   object - (xsubset0)
##   ...    - forwarded
##
## Rval: (environment)
##
environment.xsubset0 <- function (object, ...) {
  ## formula interface?
  if (is.null(object$terms)) {
    return (NULL)
  }
  
  ## extract environment
  environment(object$terms, ...)
}


## extract environment
##
## Args:
##   object - (xsubset)
##   ...    - forwarded
##
## Rval: (environment)
##
environment.xsubset <- function (object, ...) {
  NextMethod(.Generic, object, ...)
}


## extract environment
##
## Args:
##   object - (xsubset2)
##   ...    - forwarded
##
## Rval: (environment)
##
environment.xsubset2 <- function (object, ...) {
  NextMethod(.Generic, object, ...)
}


## extract variable names
##
## Args:
##   object - (xsubset)
##   size   - (integer[])
##   rank   - (integer)
##   ...    - ignored
##   full   - (logical)
##
## Rval: (character[]) variable names
##
variable.names.xsubset <- function (object, size, rank = 1, ...,
                                    full = FALSE) {
  ## trivial
  if (full) {
    return (object$x.names)
  }

  ## compute indices
  which <- object$which[, rank, size]

  ## done
  object$x.names[which]
}


## extract variable names
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - ignored
##   full   - (logical)
##
## Rval: (character[]) variable names
##
variable.names.xsubset2 <- function (object, rank = 1, ..., full = FALSE) {
  ## trivial
  if (full) {
    return (object$x.names)
  }

  ## compute indices
  which <- object$which[, rank]

  ## done
  object$x.names[which]
}


## extract formula (base method)
##
## Args:
##   x    - (xsubset0)
##   ...  - forwarded
##   full - (logical)
##
## Rval: (formula)
##
formula.xsubset0 <- function (x, ..., full = FALSE) {
  ## formula interface?
  if (is.null(x$terms)) {
    return (NULL)
  }
  
  ## extract variable names
  y.name <- x$y.name
  x.names <- variable.names(x, ..., full = full)

  ## handle intercept
  if (x$intercept) {
    x.names <- x.names[-1]
  }

  ## construct formula
  f <- paste(y.name, "~")
  f <- paste(f, paste(x.names, collapse = "+"))
  f <- paste(f, if (x$intercept) "+1" else "-1")

  ## done
  as.formula(f, environment(x, ...))
}


## extract formula
##
## Args:
##   x    - (xsubset)
##   size - (integer)
##   rank - (integer)
##   ...  - forwarded
##   full - (logical)
##
## Rval: (formula)
##
formula.xsubset <- function (x, size, rank = 1, ..., full = FALSE) {
  NextMethod(.Generic, x, size = size, rank = rank, ..., full = full)
}


## extract formula
##
## Args:
##   x    - (xsubset2)
##   rank - (integer)
##   ...  - forwarded
##   full - (logical)
##
## Rval: (formula)
##
formula.xsubset2 <- function (x, rank = 1, ..., full = FALSE) {
  NextMethod(.Generic, x, rank = rank, ..., full = full)
}


## extract model frame (base method)
##
## Args:
##   formula - (xsubset0)
##   ...     - forwarded
##   full    - (logical)
##
## Rval: (data.frame)
##
model.frame.xsubset0 <- function (formula, size, rank = 1, ..., full = FALSE) {
  ## did we keep the model?
  if (is.null(formula$model)) {
    return (NULL)
  }

  ## full model
  if (full) {
    return (formula$model)
  }

  ## extract formula
  f <- formula(formula, size = size, rank = rank, ...)

  ## done
  model.frame(f, data = formula$model, ...)
}


## extract model frame
##
## Args:
##   formula - (xsubset)
##   size    - (integer)
##   rank    - (integer)
##   ...     - forwarded
##   full    - (logical)
##
## Rval: (data.frame)
##
model.frame.xsubset <- function (formula, size, rank = 1, ..., full = FALSE) {
  NextMethod(.Generic, formula, size = size, rank = rank, ..., full = full)
}


## extract model frame
##
## Args:
##   formula - (xsubset2)
##   rank    - (integer)
##   ...     - forwarded
##   full    - (logical)
##
## Rval: (data.frame)
##
model.frame.xsubset2 <- function (formula, rank = 1, ..., full = FALSE) {
  NextMethod(.Generic, formula, rank = rank, ..., full = full)
}


## extract model matrix (base method)
##
## Args:
##   object - (xsubset0)
##   ...    - forwarded
##   full   - (logical)
##
## Rval: (matrix)
##
model.matrix.xsubset0 <- function (object, ..., full = FALSE) {
  ## did we keep the model matrix?
  if (!is.null(object[["x"]])) {
    ## full model
    if (full) {
      return (object$x)
    }

    ## extract names
    x.names <- variable.names(object, ...)

    ## done
    return (object$x[, x.names])
  }

  ## formula interface?
  if (is.null(object$terms)) {
    return (NULL)
  }
  
  ## exract formula
  f <- formula(object, ..., full = full)

  ## done
  model.matrix(f, data = model.frame(object, full = TRUE), ...)
}


## extract model matrix
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - forwarded
##   full   - (logical)
##
## Rval: (matrix)
##
model.matrix.xsubset <- function (object, size, rank = 1, ..., full = FALSE) {
  NextMethod(.Generic, object, size = size, rank = rank, ..., full = full)
}


## extract model matrix
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - forwarded
##   full   - (logical)
##
## Rval: (matrix)
##
model.matrix.xsubset2 <- function (object, rank = 1, ..., full = FALSE) {
  NextMethod(.Generic, object, rank = rank, ..., full = full)
}


## extract response variable (base method)
##
## Args:
##   object - (xsubset0)
##   ...    - forwarded
##
## Rval: (numeric[])
##
model.response.xsubset0 <- function (object, ...) {
  ## did we keep the response variable?
  if (!is.null(object[["y"]])) {
    return (object$y)
  }

  ## did we keep the model frame?
  if (is.null(object$model)) {
    return (NULL)
  }

  ## done
  model.response(object$model, ...)
}


## extract response variable
##
## Args:
##   object - (xsubset)
##   ...    - forwarded
##
## Rval: (numeric[])
##
model.response.xsubset <- function (object, ...) {
  NextMethod(.Generic, object, ...)
}


## extract response variable
##
## Args:
##   object - (xsubset2)
##   ...    - forwarded
##
## Rval: (numeric[])
##
model.response.xsubset2 <- function (object, ...) {
  NextMethod(.Generic, object, ...)
}


## extract model weights
##
## Args:
##   object - (xsubset)
##   ...    - ignored
##
## Rval: (numeric[])
##
model.weights.xsubset <- function (object, ...) {
  object$weights
}


## extract model weights
##
## Args:
##   object - (xsubset2)
##   ...    - ignored
##
## Rval: (numeric[])
##
model.weights.xsubset2 <- function (object, ...) {
  object$weights
}


## extract model offset
##
## Args:
##   object - (xsubset)
##   ...    - ignored
##
## Rval: (numeric[])
##
model.offset.xsubset <- function (object, ...) {
  object$offset
}


## extract model offset
##
## Args:
##   object - (xsubset2)
##   ...    - ignored
##
## Rval: (numeric[])
##
model.offset.xsubset2 <- function (object, ...) {
  object$offset
}


## extract terms object (base method)
##
## Args:
##   x    - (xsubset0)
##   ...  - forwarded
##   full - (logical)
##
## Rval: (terms)
##
terms.xsubset0 <- function (x, ..., full = FALSE) {
  ## formula interface?
  if (is.null(x$terms)) {
    return (NULL)
  }

  ## trivial
  if (full) {
    return (x$terms)
  }

  ## extract formula
  f <- formula(x, ...)

  ## done
  terms(f, ...)
}


## extract terms object
##
## Args:
##   x    - (xsubset)
##   size - (integer)
##   rank - (integer)
##   ...  - forwarded
##   full - (logical)
##
## Rval: (terms)
##
terms.xsubset <- function (x, size, rank = 1, ..., full = FALSE) {
  NextMethod(.Generic, x, size = size, rank = rank, ..., full = full)
}


## extract terms object
##
## Args:
##   x    - (xsubset2)
##   rank - (integer)
##   ...  - forwarded
##   full - (logical)
##
## Rval: (terms)
##
terms.xsubset2 <- function (x, rank = 1, ..., full = FALSE) {
  NextMethod(.Generic, x, rank = rank, ..., full = full)
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

  if (mask.call) {
    ## extract call to object
    cl <- match.call(expand.dots = FALSE)
    obj.call <- cl$object

    ## construct call to model frame
    mf.call <- call("model.frame")
    mf.call$formula <- obj.call
    mf.call$full <- TRUE

    ## construct call to model weights
    w.call <- call("model.weights")
    w.call$x <- obj.call

    ## construct call to model offset
    o.call <- call("model.offset")
    o.call$x <- obj.call

    ## construct call to lm
    lm.call <- call("lm")
    lm.call$formula <- f
    lm.call$data <- mf.call
    lm.call$weights <- w.call
    lm.call$offset <- o.call

    ## done
    eval(lm.call, parent.frame())
  } else {
    args <- list(formula = f,
                 data = model.frame(object, full = TRUE),
                 weights = model.weights(object),
                 offset = model.offset(object))
    do.call("lm", args)
  }
}


## refit
##
## Args:
##   object    - (xsubset)
##   size      - (integer)
##   rank      - (integer)
##   ...       - forwarded
##   mask.call - (logical) used internally
##
## Rval: (lm)
##
## Note:  'mask.call'
## Reconstruct 'call' object for pretty printing.
##
refit.xsubset <- function (object, size, rank = 1, ..., mask.call = TRUE) {
  NextMethod(.Generic, object, size = size, rank = rank, ...,
             mask.call = mask.call)
}


## refit
##
## Args:
##   object    - (xsubset2)
##   rank      - (integer)
##   ...       - forwarded
##   mask.call - (logical) used internally
##
## Rval: (lm)
##
## Note:  'mask.call'
## Reconstruct 'call' object for pretty printing.
##
refit.xsubset2 <- function (object, rank = 1, ..., mask.call = TRUE) {
  NextMethod(.Generic, object, rank = rank, ..., mask.call = mask.call)
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
  coef(refit(object, size = size, rank = rank, ..., mask.call = FALSE), ...)
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
  coef(refit(object, rank = rank, ..., mask.call = FALSE), ...)
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
  vcov(refit(object, size = size, rank = rank, ..., mask.call = FALSE), ...)
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
  vcov(refit(object, rank = rank, ..., mask.call = FALSE), ...)
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
  fitted(refit(object, size = size, rank = rank, ..., mask.call = FALSE), ...)
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
  fitted(refit(object, rank = rank, ..., mask.call = FALSE), ...)
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
  residuals(refit(object, size = size, rank = rank, ..., mask.call = FALSE), ...)
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
  residuals(refit(object, rank = rank, ..., mask.call = FALSE), ...)
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
  object$value[rank, size]
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
  object$value[rank]
}


## extract log-likelihood
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - forwarded
##
## Rval: (logLik)
##
logLik.xsubset <- function (object, size, rank = 1, ...) {
  logLik(refit(object, size = size, rank = rank, ..., mask.call = FALSE), ...)
}


## extract log-likelihood
##
## Args:
##   object - (xsubset2)
##   rank   - (integer)
##   ...    - forwarded
##
## Rval: (logLik)
##
logLik.xsubset2 <- function (object, rank = 1, ...) {
  logLik(refit(object, rank = rank, ..., mask.call = FALSE), ...)
}


## compute AIC
##
## Args:
##   object - (xsubset)
##   size   - (integer[])
##   rank   - (integer)
##   ...    - forwarded
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

  ## refits
  fitted <- lapply(size, function (sz) {
    refit(object, size = sz, rank = rank, ..., mask.call = FALSE)
  })

  ## call AIC
  rval <- do.call("AIC", c(fitted, ..., k = k))
  if (length(size) > 1) {
    row.names(rval) <- size
  } else {
    names(rval) <- size
  }

  #done
  rval
}


## compute AIC
##
## Args:
##   object - (xsubset2)
##   rank   - (integer[])
##   ...    - forwarded
##   k      - (integer)
##
## Rval: (numeric|data.frame)
##
AIC.xsubset2 <- function (object, rank = 1, ..., k = 2) {
  ## refits
  fitted <- lapply(rank, function (rk) {
    refit(object, rank = rk, ..., mask.call = FALSE)
  })

  ## call AIC
  rval <- do.call("AIC", c(fitted, ..., k = k))
  if (length(size) > 1) {
    row.names(rval) <- rank
  } else {
    names(rval) <- rank
  }

  #done
  rval
}


##
## SUMMARY
##


## summary for 'xsubset' objects
##
## Args:
##   object  - (xsubset)
##   size    - (integer[])
##   penalty - (numeric|character)
##   ...     - forwarded
##
## Rval: (summary.xsubset)
##
summary.xsubset <- function (object, size = NULL, rank = 1,
                             penalty = 2, ...) {
  paste <- function (...) base::paste(..., sep = "")

  ## default size
  if (is.null(size)) {
    size <- object$size
  }

  ## penalty
  ## ...

  ## rss
  rss <- deviance(object, size, rank, ...)

  ## aic
  if (length(size) > 1) {
    aic <- AIC(object, size = size, rank = rank, ..., k = penalty)$AIC
    best.which <- match(min(aic), aic)
  } else {
    aic <- AIC(object, size = size, rank = rank, ..., k = penalty)
    best.which <- NULL
  }

  ## names
  names(rss) <- names(aic) <- size

  ## return value
  rval <- list(call       = object$call,
               include    = object$include,
               exclude    = object$exclude,
               nbest      = object$nbest,
               size       = size,
               rank       = rank,
               which      = object$which[, rank, size],
               rss        = rss,
               penalty    = penalty,
               aic        = aic,
               best.which = best.which)
  class(rval) <- "summary.xsubset"

  ## done
  rval
}


## summary for 'xsubset2' objects
##
## Args:
##   object      - (xsubset2)
##   rank        - (integer[])
##   ...         - forwarded
##
## Rval: (summary.xsubset2)
##
summary.xsubset2 <- function (object, rank = NULL, ...) {
  paste <- function (...) base::paste(..., sep = "")

  ## extract criteria and determine best model
  if (is.null(rank)) {
    rank <- 1:object$nbest
  }

  ## rss
  rss <- lapply(rank, function (rk) {
    deviance(refit(object, rank = rk, ..., mask.call = FALSE))
  })

  ## aic
  aic <- deviance(object, rank = rank, ...)

  ## names
  names(rss) <- names(aic) <- rank

  ## return value
  rval <- list(call    = object$call,
               nbest   = object$nbest,
               include = object$include,
               exclude = object$exclude,
               penalty = object$penalty,
               rank    = rank,
               which   = object$which[, rank],
               rss     = rss,
               aic     = aic)
  class(rval) <- "summary.xsubset2"

  ## done
  rval
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
  catln("Rank: ", paste(x$rank, "."), "/", x$nbest)

  ## variable table
  catln()
  catln("Selected variables (by size):")
  which.x <- ifelse(x$which, "x", "")
  rownames(which.x)[x$include] <- paste("+", rownames(which.x)[x$include])
  rownames(which.x)[x$exclude] <- paste("-", rownames(which.x)[x$exclude])
  colnames(which.x)[x$best.which] <- paste(colnames(which.x)[x$best.which], "*")
  print(which.x, quote = FALSE)

  ## fit
  catln()
  catln("Model fit:")
  fit <- rbind(x$aic, x$rss)
  rownames(fit) <- c("AIC", "RSS")
  colnames(fit)[x$best.which] <- paste(colnames(fit)[x$best.which], "*")
  print(fit, digits = digits)
  catln()
  catln("AIC: penalty = ", x$penalty)
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
  catln()
  catln("Rank: ", paste(paste(x$rank, "."), collapse = " "), "/", x$nbest)

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
  fit <- rbind(x$aic, x$rss)
  rownames(fit) <- c("AIC", "RSS")
  colnames(fit) <- paste(colnames(fit), ".")
  print(fit, digits = digits)
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
  sub <- paste("AIC (k = ", x$penalty, ")", sep = "")

  ## color
  col <- rep(col, length.out = 2)

  ## line type
  lty <- rep(lty, length.out = 2)

  ## plot rss
  plot(x$size, x$rss, type = type[1], main = main, sub = sub,
       xlab = xlab, ylab = ylab, col = col[1], lty = lty[1], ...)

  ## plot aic
  new.old <- getOption("new")
  par(new = TRUE)
  plot(x$size, x$aic, type = type[2], xlab = "", ylab = "",
       col = col[2], lty = lty[2], axes = FALSE, ...)
  par(new = new.old)

  ## legend
  if (legend) {
    legend("topright", c("RSS", "AIC"), lty = lty, col = col, bty = "n")
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

  ## plot rss
  plot(x$rank, x$rss, type = type[1], main = main, sub = sub,
       xlab = xlab, ylab = ylab, col = col[1], lty = lty[1], ...)

  ## plot aic
  new.old <- getOption("new")
  par(new = TRUE)
  plot(x$rank, x$aic, type = type[2], xlab = "", ylab = "",
       col = col[2], lty = lty[2], axes = FALSE, ...)
  par(new = new.old)

  ## legend
  if (legend) {
    legend("topright", c("RSS", "AIC"), lty = lty, col = col, bty = "n")
  }

  ## axes
  axis(4)

  ## done
  invisible(x)
}
