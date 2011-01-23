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
  xy0 <- na.action(cbind(x0, y0))

  ## forward call
  rval <- xsubset(xy0, ...)

  ## update return value
  rval$call <- call
  rval$terms <- t
  if (!is.null(weights)) rval$weights <- w
  if (model) rval$model <- mf
  if (x.return) rval$x <- x
  if (y.return) rval$y <- y
  rval$constrasts <- constrasts
  if (!is.null(offset)) rval$offset <- o
  rval$y.name <- deparse(attr(t, "variables")[[attr(t, "response") + 1]])

  ## done
  rval
}


## method for matrix objects
##
## Args:
##   object  - (matrix) model matrix
##   y       - (numeric[]) response variable
##   ...     - arguments forwarded to 'xsubset.data.frame'
##
## Rval: (list)
##   (see 'xsubset.default')
##
## Forwards the call to 'xsubset.default'.
##
xsubset.matrix <- function (object, y = NULL, ...) {
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## process x, y
  if (is.null(y)) {
    y <- object[,  NCOL(x), drop = FALSE]
    object <- object[, -NCOL(x), drop = FALSE]
  }

  ## forward call
  rval <- xsubset.default(object, y, ...)

  ## update return value
  rval$call <- call

  # done
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
##   criterion - (numeric|character)  information criterion
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
                             size = NULL, criterion = "RSS", tolerance = 0,
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

  ## process criterion
  if (is.character(criterion)) {
    criterion.name <- toupper(criterion)
    if (criterion.name == "RSS") {
      criterion <- 0
    } else if (criterion.name == "AIC") {
      criterion <- 2
    } else if (criterion.name == "BIC") {
      criterion <- log(nobs)
    }
  } else {
    criterion.name <- "AIC"
  }

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
                 criterion = as.double(criterion),
                 tolerance = as.double(tolerance),
                 pradius   = as.integer(pradius),
                 nbest     = as.integer(nbest),
                 ## out
                 value = as.double(rep(0, nvar0 * nbest)),
                 which = as.logical(rep(0, nvar0 * nbest * nvar0)),
                 nodes = integer(1))
  C_rval <- do.call(".C", c(name = "R_select", C_args))

  ## return value
  rval <- list(call           = call,
               nobs           = nobs,
               nvar           = nvar,
               x.names        = x.names,
               y.name         = NULL,
               include        = include,
               exclude        = exclude,
               intercept      = intercept,
               size           = NULL,
               criterion.name = criterion.name,
               criterion      = criterion,
               tolerance      = NULL,
               nbest          = nbest,
               value          = NULL,
               which          = NULL,
               .nodes = C_rval$nodes)

  ## extract value & subsets
  if (criterion.name == "RSS") {
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
    class(rval) <- "xsubset"
  } else {
    rval$tolerance <- tolerance[1]
    ## aic
    rval$value <- array(C_rval$value, dim = nbest, dimnames = list(1:nbest))
    ## sel
    rval$which <- array(NA, dim = c(nvar, nbest),
                        dimnames = list(x.names, 1:nbest))
    rval$which[which0, ] <- C_rval$which[1:(nvar0 * nbest)]
    ## class
    class(rval) <- "xselect"
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
                     format.criterion(x),
                     paste(x$tolerance, collapse = " "),
                     format(x$nbest)),
                     ncol = 1)
  colnames(val) <- ""
  rownames(val) <- c("Total regressors:", "Intercept:", "Include:", "Exclude:",
                     "Size:", "Criterion:", "Tolerance:", "N best:")

  print(val, quote = FALSE)
  
  invisible(x)
}


## print method for 'xselect' objects
##
## Args:
##   x   - (xselect)
##   ... - ignored
##
## Rval: (xsubset) invisible
##
print.xselect <- function (x, ...)
{
  catln <- function (...) base::cat(..., "\n", sep = "")
  
  catln()
  catln("Call:")
  catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))
  
  val <- as.matrix(c(format(x$nvar),
                     if (x$intercept) "yes" else "no",
                     paste(x$include, collapse = " "),
                     paste(x$exclude, collapse = " "),
                     format.criterion(x),
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
##   legend - (logical)
##   xlab   - (character)
##   ylab   - (character)
##   main   - (character)
##   col    - (integer[]|character[])
##   lty    - (integer)
##   type   - (character)
##   ...    - passed ot 'plot.default'
##
## Rval: (xsubset) invisible
##
## Graphical arguments are passed to 'plot.default'.
##
plot.xsubset <- function (x, size = NULL, rank = 1, legend = TRUE,
  xlab = NULL, ylab = "", main = NULL, col = "blue", lty = 1, type = "b", ...)
{
  if (is.null(size)) {
    size <- x$size
  }

  rss <- deviance(x, size = size, rank = rank)

  if (is.null(xlab)) {
    xlab <- "Number of regressors in model"
  }
  if (is.null(main)) {
    main <- paste("RSS (rank = ", rank, ")", sep = "")
  }

  plot(size, rss,
       ylab = ylab, xlab = xlab, main = main,
       type = type, lty = lty, col = col, ...)

  if (legend) {
    legend("topright", "RSS", lty = lty, col = col, bty = "n")
  }

  invisible(x)
}


## plot method for 'xselect' objects
##
## Args:
##   x      - (xselect)
##   rank   - (integer[])
##   legend - (logical)
##   xlab   - (character)
##   ylab   - (character)
##   main   - (character)
##   col    - (integer[]|character[])
##   lty    - (integer)
##   type   - (character)
##   ...    - passed ot 'plot.default'
##
## Rval: (xsubset) invisible
##
## Graphical arguments are passed to 'plot.default'.
##
plot.xselect <- function (x, rank = NULL, legend = TRUE,
  xlab = NULL, ylab = "", main = NULL, col = "blue", lty = 1, type = "b", ...)
{
  if (is.null(rank)) {
    rank <- 1:x$nbest
  }
  
  aic <- deviance(x, rank = rank)

  if (is.null(xlab)) {
    xlab <- "Rank"
  }
  if (is.null(main)) {
    main <- "AIC"
  }

  plot(rank, aic,
       ylab = ylab, xlab = xlab, main = main,
       type = type, lty = lty, col = col, ...)

  size <- apply(x$which[, rank], 2, sum)
  text(rank, aic, labels = size, adj = c(0, 1.5), cex = 0.8)

  if (legend) {
    legend("topleft", format.criterion(x), lty = lty, col = col, bty = "n")
  }

  invisible(x)
}


## extract environment
##
## Args:
##   object - (xsubset)
##   ...    - passed to 'environment.default'
##
## Rval: (environment)
##
environment.xsubset <- function (object, ...) {
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
##   object - (xselect)
##   ...    - passed to 'environment.default'
##
## Rval: (environment)
##
environment.xselect <- function (object, ...) {
  ## formula interface?
  if (is.null(object$terms)) {
    return (NULL)
  }
  
  ## extract environment
  environment(object$terms, ...)
}


## extract variable names
##
## Args:
##   object - (xsubset)
##   full   - (logical)
##   size   - (integer[])
##   rank   - (integer)
##   ...    - ignored
##
## Rval: (character[]) variable names
##
variable.names.xsubset <- function (object, full = FALSE, size = NULL,
                                    rank = 1, ...) {
  ## trivial
  if (full) {
    return (object$x.names)
  }

  ## default size
  if (is.null(size)) {
    size <- object$size
  }

  ## compute indices
  which <- object$which[, rank, size]

  ## done
  object$x.names[which]
}


## extract variable names
##
## Args:
##   object - (xselect)
##   full   - (logical)
##   rank   - (integer)
##   ...    - ignored
##
## Rval: (character[]) variable names
##
variable.names.xselect <- function (object, full = FALSE, rank = 1, ...) {
  ## trivial
  if (full) {
    return (object$x.names)
  }

  ## compute indices
  which <- object$which[, rank]

  ## done
  object$x.names[which]
}


## extract formula
##
## Args:
##   x    - (xsubset)
##   full - (logical)
##   size - (integer)
##   rank - (integer)
##   ...  - passed to 'variable.names.xsubset'
##
## Rval: (formula)
##
formula.xsubset <- function (x, full = FALSE, size = NULL, rank = 1, ...) {
  ## formula interface?
  if (is.null(x$terms)) {
    return (NULL)
  }
  
  ## extract variable names
  y.name <- x$y.name
  x.names <- variable.names(x, full = full, size = size, rank = rank, ...)

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
##   x    - (xselect)
##   full - (logical)
##   rank - (integer)
##   ...  - passed to 'variable.names.xselect'
##
## Rval: (formula)
##
formula.xselect <- function (x, full = FALSE, rank = 1, ...) {
  ## formula interface?
  if (is.null(x$terms)) {
    return (NULL)
  }
  
  ## extract variable names
  y.name <- x$y.name
  x.names <- variable.names(x, full = full, rank = rank, ...)

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


## extract model frame
##
## Args:
##   formula - (xsubset)
##   full    - (logical)
##   size    - (integer)
##   rank    - (integer)
##   ...     - passed to 'formula.xsubset'
##
## Rval: (data.frame)
##
model.frame.xsubset <- function (formula, full = FALSE, size = NULL, rank = 1, ...) {
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
  model.frame(f, data = formula$model)
}


## extract model frame
##
## Args:
##   formula - (xselect)
##   full    - (logical)
##   rank    - (integer)
##   ...     - passed to 'formula.xselect'
##
## Rval: (data.frame)
##
model.frame.xselect <- function (formula, full = FALSE, rank = 1, ...) {
  ## did we keep the model?
  if (is.null(formula$model)) {
    return (NULL)
  }

  ## full model
  if (full) {
    return (formula$model)
  }

  ## extract formula
  f <- formula(formula, rank = rank, ...)

  ## done
  model.frame(f, data = formula$model)
}


## extract model matrix
##
## Args:
##   object - (xsubset)
##   full   - (logical)
##   size   - (integer)
##   rank   - (integer)
##   ...    - ignored
##
## Rval: (matrix)
##
model.matrix.xsubset <- function (object, full = FALSE, size = NULL, rank = 1, ...) {
  ## did we keep the model matrix?
  if (!is.null(object$x)) {
    ## full model
    if (full) {
      return (object$x)
    }

    ## extract names
    x.names <- variable.names(object, size = size, rank = rank)

    ## done
    return (object$x[, x.names])
  }

  ## formula interface?
  if (is.null(object$terms)) {
    return (NULL)
  }
  
  ## exract formula
  f <- formula(object, size = size, rank = rank)

  ## done
  model.matrix(f, data = model.frame(object, full = TRUE))
}


## extract model matrix
##
## Args:
##   object - (xselect)
##   full   - (logical)
##   rank   - (integer)
##   ...    - ignored
##
## Rval: (matrix)
##
model.matrix.xselect <- function (object, full = FALSE, rank = 1, ...) {
  ## did we keep the model matrix?
  if (!is.null(object$x)) {
    ## full model
    if (full) {
      return (object$x)
    }

    ## extract names
    x.names <- variable.names(object, rank = rank)

    ## done
    return (object$x[, x.names])
  }

  ## formula interface?
  if (is.null(object$terms)) {
    return (NULL)
  }

  ## exract formula
  f <- formula(object, rank = rank)

  ## done
  model.matrix(f, data = model.frame(object, full = TRUE))
}


## extract response variable
##
## Args:
##   object - (xsubset)
##   ...    - passed to 'model.response'
##
## Rval: (numeric[])
##
model.response.xsubset <- function (object, ...) {
  ## did we keep the response variable?
  if (!is.null(object$y)) {
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
##   object - (xselect)
##   ...    - passed to 'model.response'
##
## Rval: (numeric[])
##
model.response.xselect <- function (object, ...) {
  ## did we keep the response variable?
  if (!is.null(object$y)) {
    return (object$y)
  }

  ## did we keep the model frame?
  if (is.null(object$model)) {
    return (NULL)
  }

  ## done
  model.response(object$model, ...)
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
##   object - (xselect)
##   ...    - ignored
##
## Rval: (numeric[])
##
model.weights.xselect <- function (object, ...) {
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
##   object - (xselect)
##   ...    - ignored
##
## Rval: (numeric[])
##
model.offset.xselect <- function (object, ...) {
  object$offset
}


## extract terms object
##
## Args:
##   x    - (xsubset)
##   full - (logical)
##   size - (integer)
##   rank - (integer)
##   ...  - passed to 'terms.formula'
##
## Rval: (terms)
##
terms.xsubset <- function (x, full = FALSE, size = NULL, rank = 1, ...) {
  ## formula interface?
  if (is.null(x$terms)) {
    return (NULL)
  }

  ## trivial
  if (full) {
    return (x$terms)
  }

  ## extract formula
  f <- formula(x, size = size, rank = rank)

  ## done
  terms(f, ...)
}


## extract terms object
##
## Args:
##   x    - (xselect)
##   full - (logical)
##   size - (integer)
##   rank - (integer)
##   ...  - passed to 'terms.formula'
##
## Rval: (terms)
##
terms.xselect <- function (x, full = FALSE, rank = 1, ...) {
  ## formula interface?
  if (is.null(x$terms)) {
    return (NULL)
  }

  ## trivial
  if (full) {
    return (x$terms)
  }

  ## extract formula
  f <- formula(x, rank = rank)

  ## done
  terms(f, ...)
}


## refit associated 'lm' object
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - ignored
##   mask.call - (logical) used internally
##
## Rval: (lm)
##
## Note:  'mask.call'
## Reconstruct 'call' object for pretty printing.
##
refit.xsubset <- function (object, size = NULL, rank = 1, ..., mask.call = TRUE) {
  ## formula interface?
  if (is.null(object$terms)) {
    return (NULL)
  }

  ## extract formula
  f <- formula(object, size = size, rank = rank)

  if (mask.call) {
    ## extract call to object
    cl <- match.call(expand.dots = FALSE)
    obj.call <- cl$object

    ## construct call to model frame
    mf.call <- call("model.frame")
    mf.call$formula <- obj.call

    ## construct call to model weights
    w.call <- call("model.weights")
    w.call$object <- obj.call

    ## construct call to model offset
    o.call <- call("model.offset")
    o.call$object <- obj.call

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


## refit associated 'lm' object
##
## Args:
##   object - (xselect)
##   size   - (integer)
##   rank   - (integer)
##   ...    - ignored
##   mask.call - (logical) used internally
##
## Rval: (lm)
##
## Note:  'mask.call'
## Reconstruct 'call' object for pretty printing.
##
refit.xselect <- function (object, rank = 1, ..., mask.call = TRUE) {
  ## formula interface?
  if (is.null(object$terms)) {
    return (NULL)
  }

  ## extract formula
  f <- formula(object, rank = rank)

  if (mask.call) {
    ## extract call to object
    cl <- match.call(expand.dots = FALSE)
    obj.call <- cl$object

    ## construct call to model frame
    mf.call <- call("model.frame")
    mf.call$formula <- obj.call

    ## construct call to model weights
    w.call <- call("model.weights")
    w.call$object <- obj.call

    ## construct call to model offset
    o.call <- call("model.offset")
    o.call$object <- obj.call

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


## extract ceofficients
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - passed to 'coef.lm'
##
## Rval: (numeric[])
##
coef.xsubset <- function (object, size = NULL, rank = 1, ...) {
  coef(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract ceofficients
##
## Args:
##   object - (xselect)
##   size   - (integer)
##   rank   - (integer)
##   ...    - passed to 'coef.lm'
##
## Rval: (numeric[])
##
coef.xselect <- function (object, rank = 1, ...) {
  coef(refit(object, rank = rank, mask.call = FALSE), ...)
}


## extract variance-covariance matrix
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - passed to 'vcov.lm'
##
## Rval: (matrix)
##
vcov.xsubset <- function (object, size = NULL, rank = 1, ...) {
  vcov(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract variance-covariance matrix
##
## Args:
##   object - (xselect)
##   size   - (integer)
##   rank   - (integer)
##   ...    - passed to 'vcov.lm'
##
## Rval: (matrix)
##
vcov.xselect <- function (object, rank = 1, ...) {
  vcov(refit(object, rank = rank, mask.call = FALSE), ...)
}


## extract fitted values
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - passed to 'fitted.lm'
##
## Rval: (numeric[])
##
fitted.xsubset <- function (object, size = NULL, rank = 1, ...) {
  fitted(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract fitted values
##
## Args:
##   object - (xselect)
##   size   - (integer)
##   rank   - (integer)
##   ...    - passed to 'fitted.lm'
##
## Rval: (numeric[])
##
fitted.xselect <- function (object, rank = 1, ...) {
  fitted(refit(object, rank = rank, mask.call = FALSE), ...)
}


## extract residuals
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - passed to 'residuals.lm'
##
## Rval: (numeric[])
##
residuals.xsubset <- function (object, size = NULL, rank = 1, ...) {
  residuals(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract residuals
##
## Args:
##   object - (xselect)
##   size   - (integer)
##   rank   - (integer)
##   ...    - passed to 'residuals.lm'
##
## Rval: (numeric[])
##
residuals.xselect <- function (object, rank = 1, ...) {
  residuals(refit(object, rank = rank, mask.call = FALSE), ...)
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
## Returns the deviance for the specified
## subset size and rank.
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
## Returns the deviance for the specified
## subset size and rank.
##
deviance.xselect <- function (object, rank = 1, ...) {
  object$value[rank]
}


## extract log-likelihood
##
## Args:
##   object - (xsubset)
##   size   - (integer)
##   rank   - (integer)
##   ...    - passed to 'logLik.lm'
##
## Rval: (logLik)
##
logLik.xsubset <- function (object, size = NULL, rank = 1, ...) {
  logLik(refit(object, size = size, rank = rank, mask.call = FALSE), ...)
}


## extract log-likelihood
##
## Args:
##   object - (xselect)
##   size   - (integer)
##   rank   - (integer)
##   ...    - passed to 'logLik.lm'
##
## Rval: (logLik)
##
logLik.xselect <- function (object, rank = 1, ...) {
  logLik(refit(object, rank = rank, mask.call = FALSE), ...)
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
    refit(object, size = sz, rank = rank, mask.call = FALSE)
  })

  ## call AIC
  rval <- do.call("AIC", c(fitted, k = k))
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
##   object - (xselect)
##   rank   - (integer[])
##   ...    - ignored
##   k      - (integer)
##
## Rval: (numeric|data.frame)
##
## Note:  'size'
## The default behavior is to extract the AIC for
## all selected subset sizes.
##
AIC.xsubset <- function (object, rank = 1, ..., k = 2) {
  ## refits
  fitted <- lapply(rank, function (rk) {
    refit(object, rank = rk, mask.call = FALSE)
  })

  ## call AIC
  rval <- do.call("AIC", c(fitted, k = k))
  if (length(size) > 1) {
    row.names(rval) <- size
  } else {
    names(rval) <- size
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
##   object      - (xsubset)
##   size        - (integer[])
##   criterion   - (numeric|character)
##   ...         - ignored
##
## Rval: (summary.xsubset)
##
summary.xsubset <- function (object, size = NULL, rank = 1,
                             criterion = "AIC", ...) {
  paste <- function (...) base::paste(..., sep = "")

  ## default size
  if (is.null(size)) {
    size <- object$size
  }

  ## criterion
  if (is.character(criterion)) {
    criterion.name <- toupper(criterion)
    if (criterion.name == "AIC") {
      criterion <- 2
    } else if (criterion.name == "BIC") {
      criterion <- log(object$nobs)
    }
  } else {
    criterion.name <- "AIC"
  }

  ## rss
  rss <- deviance(object, size, rank)

  ## construct variable table
  which <- matrix(FALSE, nrow = object$nvar, ncol = length(size))
  dimnames(which) <- list(variable.names(object, full = TRUE), size)
  rownames(which)[object$include] <- paste("+", rownames(which)[object$include])
  rownames(which)[object$exclude] <- paste("-", rownames(which)[object$exclude])
  colnames(which)[best.which] <- paste(size[best.which], "*")
  which[, size] <- object$which[, rank, size]

  ## aic
  if (length(size) > 1) {
    aic <- AIC(object, size = size, rank = rank, k = criterion)$AIC
    best.which <- match(min(aic), aic)
  } else {
    aic <- AIC(object, size = size, rank = rank, k = criterion)
    best.which <- NULL
  }

  ## names
  names(rss) <- names(aic) <- size

  ## return value
  rval <- list(call           = object$call,
               include        = object$include,
               exclude        = object$exclude,
               nbest          = object$nbest,
               size           = size,
               rank           = rank,
               which          = object$which[, rank, size],
               rss            = rss,
               criterion.name = criterion.name,
               criterion      = criterion,
               aic            = aic,
               best.which     = best.which)
  class(rval) <- "summary.xsubset"

  ## done
  rval
}


## summary for 'xselect' objects
##
## Args:
##   object      - (xselect)
##   rank        - (integer[])
##   ...         - ignored
##
## Rval: (summary.xselect)
##
summary.xselect <- function (object, rank = NULL, ...) {
  paste <- function (...) base::paste(..., sep = "")

  ## extract criteria and determine best model
  if (is.null(rank)) {
    rank <- 1:object$nbest
  }
  criterion <- object$criterion.raw

  ## rss
  rss <- lapply(rank, function (rk) {
    deviance(refit(object, rank = rk, mask.call = FALSE))
  })

  ## aic
  aic <- deviance(object, rank = rank)

  ## names
  names(rss) <- names(aic) <- rank

  ## return value
  rval <- list(call           = object$call,
               nbest          = object$nbest,
               include        = object$include,
               exclude        = objcet$exclude,
               criterion.name = object$criterion.name,
               criterion      = object$criterion,
               rank           = rank,
               which          = object$which[, rank],
               rss            = rss,
               aic            = aic)
  class(rval) <- "summary.xselect"

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
  catln <- function (...) base::cat(..., "\n", sep = "")

  ## call
  catln()
  catln("Call:")
  catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))

  ## rank
  catln()
  catln("Rank: ", paste(x$rank, ".", sep = ""), "/", x$nbest)

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
  rownames(fit) <- c(pretty.criterion(x), "RSS")
  print(fit, digits = digits)
  catln()
  catln("* best subset")
    
  invisible(x)
}


## print 'xselect' summary
##
## Args:
##   x      - (summary.xselect)
##   digits - (integer)
##   ...    - ignored
##
## Rval: (summary.xselect) invisible
##
print.summary.xselect <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
  catln <- function (...) base::cat(..., "\n", sep = "")

  ## call
  catln()
  catln("Call:")
  catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))

  ## rank
  catln()
  catln("Rank: ", paste(x$rank, ".", sep = ""), "/", x$nbest)

  ## variables
  catln()
  catln("Selected variables (by rank):")
  which.x <- ifelse(x$which, "x", "")
  rownames(which.x)[x$include] <- paste("+", rownames(which.x)[x$include])
  rownames(which.x)[x$exclude] <- paste("-", rownames(which.x)[x$exclude])
  print(which.x, quote = FALSE)

  ## fit
  catln()
  catln("Model fit:")
  fit <- rbind(x$aic, x$rss)
  rownames(fit) <- c(pretty.criterion(x), "RSS")
  print(fit, digits = digits)
    
  ## done
  invisible(x)
}


## plot 'xsubset' summary
##
## Args:
##   x      - (xsubset)
##   size   - (integer)
##   legend - (logical)
##   xlab   - (character)
##   ylab   - (character)
##   main   - (character)
##   col    - (integer[]|character[])
##   lty    - (integer)
##   type   - (character)
##   ...    -
##
## Rval: (summary.xsubset) invisible
##
## All arguments are passed to 'plot.default'.
##
plot.summary.xsubset <- function (x, legend = TRUE, xlab = NULL, ylab = "",
  main = NULL, col = c("blue", "red"), lty = 1, type = "b", ...)
{
  ## x label
  if (is.null(xlab)) {
    xlab <- "Number of regressors in model"
  }

  ## main title
  if (is.null(main)) {
    main <- paste("AIC and residual sum of squares")
  }

  ## color
  col <- rep(col, length.out = 2)

  ## type
  type <- rep(type, length.out = 2)

  ## line type
  lty <- rep(lty, length.out = 2)

  ## plot rss
  plot(x$size, x$rss, xlab = xlab, ylab = ylab, main = main,
    col = col[1], lty = lty[1], type = type[1], ...)

  ## plot aic
  new.old <- getOption("new")
  par(new = TRUE)
  plot(x$size, x$aic, xlab = "", ylab = "", col = col[2],
       lty = lty[2], type = type[2], axes = FALSE, ...)
  par(new = new.old)

  ## legend
  if (legend) {
    legend("topright", c("RSS", pretty.criterion(x)),
           lty = lty, col = col, bty = "n")
  }

  ## axes
  axis(4)

  ## done
  invisible(x)
}


## plot 'xselect' summary
##
## Args:
##   x      - (xselect)
##   legend - (logical) legend?
##   xlab   - (character) x label
##   ylab   - (character) y label
##   main   - (character) main title
##   col    - (integer[]|character[]) color
##   lty    - (integer) line type
##   type   - (character) plot type
##   ...    -
##
## Rval: (summary.xselect) invisible
##
## All arguments are passed to 'plot.default'.
##
plot.summary.xselect <- function (x, legend = TRUE, xlab = NULL, ylab = "",
  main = NULL, col = c("blue", "red"), lty = 1, type = "b", ...)
{
  ## x label
  if (is.null(xlab)) {
    xlab <- "Number of regressors in model"
  }

  ## main title
  if (is.null(main)) {
    main <- paste("AIC and residual sum of squares")
  }

  ## color
  col <- rep(col, length.out = 2)

  ## type
  type <- rep(type, length.out = 2)

  ## line type
  lty <- rep(lty, length.out = 2)

  ## plot rss
  plot(x$rank, x$rss, xlab = xlab, ylab = ylab, main = main,
    col = col[1], lty = lty[1], type = type[1], ...)

  ## plot aic
  new.old <- getOption("new")
  par(new = TRUE)
  plot(x$rank, x$aic, xlab = "", ylab = "", col = col[2],
       lty = lty[2], type = type[2], axes = FALSE, ...)
  par(new = new.old)

  ## legend
  if (legend) {
    legend("topright", c("RSS", pretty.criterion(x)),
           lty = lty, col = col, bty = "n")
  }

  ## axes
  axis(4)

  ## done
  invisible(x)
}



##
## MISC
##


format.criterion <- function (x) {
  if (x$criterion.name %in% c("RSS", "BIC")) {
    x$criterion.name
  } else {
    paste("AIC (k = ", format(x$criterion), ")", sep = "")
  }
}
