##
## File:  xsubset.R
##
## Author:  Achim Zeileis
##
## Modified by:  Marc Hofmann
##


## method for matrix objects
##
## Args:
##   object  - (matrix) model matrix
##   y       - (~matrix = NULL) response variable
##   ... - arguments forwarded to 'xsubset.data.frame'
##
## Rval: (list)
##   (see 'xsubset.default')
##
## This method coerces the matrix 'cbind(y, x)' to a
## 'data.fame' and forwards the call to 'xsubset.data.frame'.
##
## Note:  response variable
## If 'y' is 'NULL', the last column of 'x' is used as the
## response variable.
##
## Note:  intercept
## No automatic intercept detection.
## 
##
xsubset.matrix <- function (object, y = NULL, ...) {
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## process x, y
  x <- object
  if (is.null(y)) {
    y <- x[, ncol(x), drop = FALSE]
    x <- x[, -ncol(x), drop = FALSE]
  }
  y <- as.matrix(y)
  if (!is.numeric(x)) {
    stop ("'object' (model matrix) must be numeric")
  }
  if (!is.numeric(y)) {
    stop ("'y' must be numeric")
  }
  if (nrow(x) != nrow(y)) {
    stop ("'object' (model matrix) and 'y' are non-conforming")
  }
  if (ncol(y) > 1) {
    warning ("dropping redundant response variables")
    y <- y[, 1, drop = FALSE]
  }
  nobs <- nrow(x)
  nreg <- ncol(x)

  ## construct data frame
  df <- as.data.frame(cbind(y, x))

  ## forward call
  rval <- xsubset(df, ...)

  ## update return value
  rval$call <- call

  # done
  rval
}


## method for data fame objects
##
## Args:
##   object  - (data.frame) data frame
##   y - (integer|character = 1) index of response variable
##   ... - arguments forwarded to 'xsubset.formula'
##
## Rval: (list)
##   (see 'xsubset.default')
##
## This method simply constructs a 'formula' from the
## object ('formula(object, env = NULL)') and forwards
## the call to 'xsubset.formula'.
##
xsubset.data.frame <- function (object, y = 1, ...) {
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## rearrange data frame
  if (is.character(y)) {
    y <- match(y, names(object))
  }
  df <- data.frame(object[y], object[-y])

  ## extract formula
  f <- formula(df, env = NULL)

  ## forward call
  rval <- xsubset(f, object, ...)

  ## update return value
  rval$call <- call

  # done
  rval
}


## standard formula interface
##
## Args:
##   object - (formula) 'lm''s 'formula' argument
##   data - (data.frame = NULL)
##   row.subset - (vector = NULL) 'lm''s 'subset' argument
##   weights - (numeric[] = NULL)
##   na.action - (function)
##   model - (logical)
##   x - (logical)
##   y - (logical)
##   contrasts - (numeric[])
##   offset - (numeric[])
##   ... - arguments forwarded to 'xsubset.lm'
##
## Rval: (list)
##   (see 'xsubset.default')
##
## Constructs an 'lm' object and forwards the call
## to 'xsubset.lm'.
##
## Note:  arguments
## The named arguments are passed to 'lm'.
##
xsubset.formula <- function (object, data = NULL, row.subset = NULL, weights = NULL,
                             na.action = na.omit, model = TRUE, x = FALSE, y = FALSE,
                             contrasts = NULL, offset = NULL, ...)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## construct lm object
  ##
  ## TODO:  'keep.order'
  ## Is it necessary to freeze the order of the terms
  ## by way of 'terms(formula, keep.order = TRUE)'?
  ##
  cl <- match.call(expand.dots = FALSE)
  ma <- match(c("object", "data", "row.subset", "weights", "na.action",
                "model", "x", "y", "contrasts", "offset"), names(cl), 0)
  cl <- cl[c(1, ma)]
  cl$formula <- cl$object
  cl$object <- NULL
  cl$subset <- cl$row.subset
  cl$row.subset <- NULL
  cl[[1]] <- as.name("lm")
  lm <- eval(cl, parent.frame())

  ## forward call
  rval <- xsubset(lm, ...)

  ## update return value
  rval$call <- call

  ## done
  rval
}


## interface for fitted lm regression
##
## Args:
##   lm - (lm)
##   ... - arguments forwarded to 'xsubset.default'
##
## Rval:
##   (see 'xsubset.default')
##
## Does nothing.  Simply forwards to 'xsubset.default'.
##   
xsubset.lm <- function (object, ...)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## forward call
  rval <- xsubset.default(object, ...)

  ## update return value
  rval$call <- call
  
  ## done
  rval
}


## workhorse method
##
## Args:
##   object - (lm)
##   include - (integer[] = NULL)
##   exclude - (integer[] = NULL)
##   size - (integer[] = NULL)
##   tolerance - (numeric[] = 0)
##   pradius - (integer)
##   precision - (character)
##   ... - ignored
##
## Rval: (list)
##   call - (call)
##   lm   - (lm)
##   nreg - (integer)
##   include - (integer[])
##   exclude - (integer[])
##   size - (integer[])
##   which - (list)
##   rss - (numeric[])
##   intercept - (logical)
##
## All the fun happens here.
##
xsubset.default <- function (object, include = NULL, exclude = NULL,
                             size = NULL, tolerance = 0, ...,
                             pradius = NULL, precision = "double")
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## extract information from lm object
  x <- model.matrix(object)
  y <- model.response(object)

  ## dimensions
  nobs <- nrow(x)
  nreg <- ncol(x)

  ## extract weights and offset
  w <- model.weights(object)
  o <- model.offset(object)

  ## handle weights
  if (is.null(w)) {
    w <- rep(1, NROW(x))
  }
  wnz <- w != 0
  x.computed <- sqrt(w[wnz]) * x[wnz, , drop = FALSE]
  y.computed <- sqrt(w[wnz]) * y[wnz, drop = FALSE]

  ## handle offset
  if (is.null(o)) {
    o <- rep(0, nobs)
  }
  y.computed <- y.computed - o[wnz]

  ## dimensions
  nobs.computed <- nrow(x.computed)

  ## variables names
  x.names <- variable.names(object)

  ## include processing
  include.computed <- include
  if (is.null(include.computed)) {
    include.computed <- numeric(0)
  }
  ## character --> numeric
  if (is.character(include.computed)) {
    include.computed <- match(include.computed, x.names)
  }
  ## logical --> numeric
  if (is.logical(include.computed)) {
    include.computed <- which(rep(include.computed, length.out = nreg))
  }
  ## remove NAs
  if (any(is.na(include.computed))) {
    warning ("omitting non-existing columns in 'include'")
    include.computed <- na.omit(include.computed)
  }
  ## remove non-positive indices
  if (any(include.computed <= 0)) {
    warning ("omitting non-positive indexes in 'include'")
    include.computed <- include.computed[include.computed > 0]
  }
  ## handle index overflow
  if (any(include.computed > nreg)) {
    warning ("'include' indexes too large; fixing 'include'")
    include.computed <- include.computed[include.computed <= nreg]
  }
  ## canonicalize include
  include.computed <- sort(unique(round(include.computed)))

  ## exclude processing
  exclude.computed <- exclude
  if (is.null(exclude.computed)) {
    exclude.computed <- numeric(0)
  } else {
    ## character --> numeric
    if (is.character(exclude.computed)) {
      exclude.computed <- match(exclude.computed, x.names)
    }
    ## logical --> numeric
    if (is.logical(exclude.computed)) {
      exclude.computed <- which(rep(exclude.computed, length.out = nreg))
    }
    ## remove NAs
    if (any(is.na(exclude.computed))) {
      warning ("non-existing columns selected in 'exclude'")
      exclude.computed <- na.omit(exclude.computed)
    }
    ## handle index underflow
    if (any(exclude.computed <= 0)) {
      warning ("'exclude' indexes must be positive")
      exclude.computed <- exclude.computed[exclude.computed > 0]
    }  
    ## handle index overflow
    if (any(exclude.computed > nreg)) {
      warning ("'exclude' indexes too large; fixing 'exclude'")
      exclude.computed <- exclude.computed[exclude.computed <= nreg]
    }
    ## canonicalize indices
    exclude.computed <- sort(unique(round(exclude.computed)))
  }

  ## include/exclude non-overlapping
  if (any(intersect(include.computed, exclude.computed))) {
    warning ("'include' and 'exclude' must not overlap; fixing 'exclude'")
    exclude.computed <- setdiff(exclude.computed, include.computed)
  }

  ## handle intercept
  intercept.computed <- x.names[1] == "(Intercept)"
  if (intercept.computed) {
    if (is.element(1, include.computed)) {
      ## OK, already selected
    } else if (is.element(1, exclude.computed)) {
      ## not selected
      intercept.computed <- FALSE
    } else {
      ## select
      include.computed <- c(1, include.computed)
    }
  }

  ## include/exclude columns
  free <- setdiff(1:nreg, c(include.computed, exclude.computed))
  column.index <- c(include.computed, free)
  x.computed <- x.computed[, column.index, drop = FALSE]
  nreg.computed <- ncol(x.computed)

  ## mark included columns      
  mark.computed <- length(include.computed)

  ## process size
  size.computed <- process.size(size, size.default = (mark.computed + 1):nreg.computed,
                                size.allow = (mark.computed + 1):nreg.computed,
                                allow.vector = TRUE, canonicalize = TRUE)
  nsize <- length(size.computed)
  size.mask <- 1:nreg.computed %in% size.computed

  ## process tolerance
  if (!is.numeric(tolerance)) {
    stop ("'tolerance' must be numeric")
  }
  if (length(tolerance) > nsize) {
    warning ("redundant values in 'tolerance'; truncating 'tolerance'")
    tolerance <- tolerance[1:nsize]
  }
  tolerance <- rep(tolerance, length.out = nsize)
  tolerance.computed <- rep(0, nreg.computed)
  tolerance.computed[size.mask] <- tolerance
  tolerance.computed[!size.mask] <- .Machine$double.xmax

  ## process pradius
  pradius.computed <- pradius
  if (is.null(pradius.computed)) {
    pradius.computed <- round(nreg.computed/3)
  }
  if (!is.numeric(pradius.computed)) {
    stop ("'pradius' must be numeric")
  }
  if (length(pradius.computed) > 1) {
    pradius <- pradius.computed[1]
    warning ("redundant values in 'pradius'; using ", pradius.computed)
  }
  pradius.rounded <- round(pradius.computed)
  if (pradius.computed != pradius.rounded) {
    warning ("'pradius' must be integer; rounding 'pradius'")
    pradius.computed <- pradius.rounded
  }

  ## precision
  precision.computed <- precision
  if (length(precision.computed) > 1) {
    precision.computed <- precision.computed[1]
    warning ("redundant values in 'precision'; using ", precision.computed)
  }
  if (!(precision.computed %in% c("single", "double"))) {
    stop ("'precision' must be \"single\" or \"double\"")
  }

  ## call underlying C code
  C_args <- list(## in
                 precision = as.character(precision.computed),
                 nobs      = as.integer(nobs.computed),
                 nreg      = as.integer(nreg.computed),
                 xy        = as.double(cbind(x.computed, y.computed)),
                 mark      = as.integer(mark.computed),
                 tolerance = as.double(tolerance.computed),
                 pradius   = as.integer(pradius.computed),
                 ## out
                 rss   = double(nreg.computed),
                 which = integer(nreg.computed * nreg.computed),
                 count = integer(1))
  C_rval <- do.call(".C", c(name = "R_subset", C_args))


  ## extract selected regressors
  which <- lapply(size.computed, function (sz) {
    first <- (sz - 1) * nreg.computed + 1
    last <- first + sz - 1
    which <- C_rval$which[first:last] + 1
    column.index[which]
  })

  ## extract RSS
  rss <- C_rval$rss[size.computed]
  
  ## names
  names(which) <- names(rss) <- size.computed

  # return value
  rval <- list(call      = call,
               lm        = object,
               nreg      = nreg,
               include   = include.computed,
               exclude   = exclude.computed,
               size      = size.computed,
               which     = which,
               rss       = rss,
               intercept = intercept.computed,
               .nodes = C_rval$count)
  class(rval) <- "xsubset"

  # done
  rval
}


## print method for 'xsubset' objects
##
## Args:
##   x - (xsubset)
##
## Rval: (xsubset) invisible
##
print.xsubset <- function (x, ...)
{
  catln <- function (...) base::cat(..., "\n", sep = "")
  
  catln()
  catln("Call:")
  catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))

  ## format include, exclude, size
  size.pretty <- pretty.integer(x$size)
  include.pretty <- pretty.integer(x$include)
  exclude.pretty <- pretty.integer(x$exclude)

  ## format nreg and intercept
  nreg.pretty <- format(x$nreg)
  intercept.pretty <- if (x$intercept) "YES" else "NO"
  
  val <- as.matrix(c(nreg.pretty, include.pretty, exclude.pretty,
                     size.pretty, intercept.pretty), ncol = 1)
  colnames(val) <- ""
  rownames(val) <- c("Number of regressors:","Include:", "Exclude:",
                     "Subset sizes assessed:", "Intercept:")
  print(val, quote = FALSE)
  
  invisible(x)
}


## plot method for 'xsubset' objects
##
## Args:
##   x - (xsubset)
##   size - (integer = NULL)
##   legend - (logical = TRUE)
##   xlab - (character = ...)
##   ylab - (character = "")
##   main - (character = NULL)
##   col - (integer|character = "blue")
##   lty - (integer = 1)
##   type - (character = "b")
##   ... - passed ot 'plot.default'
##
## Rval: (xsubset) invisible
##
## All arguments are passed to 'plot.default'.
##
plot.xsubset <- function (x, size = NULL, legend = TRUE,
  xlab = "Number of regressors in model", ylab = "",
  main = NULL, col = "blue", lty = 1, type = "b", ...)
{
  ## process size
  size <- process.size(size, size.default = x$size, size.allow = x$size,
                       allow.vector = TRUE)

  ## extract rss
  rss <- deviance(x, size = size)

  ## title
  if (is.null(main)) {
    main <- paste("Residual sum of squares")
  }

  ## plot
  plot(size, as.vector(rss),
       ylab = ylab, xlab = xlab, main = main,
       type = type, lty = lty, col = col, ...)

  ## legend
  if (legend) {
    legend("topright", "RSS", lty = lty, col = col, bty = "n")
  }

  ## axis
  axis(4)
  
  ## done
  invisible(x)
}


## extract environment
##
## Args:
##   object - (xsubset)
##   ... - passed to 'environment.default'
##
## Rval: (environment)
##
## Returns the 'environment' object associated with
## this object's 'lm' object.
##
environment.xsubset <- function (object, ...) {
  environment(object$lm, ...)
}


## extract variable names
##
## Args:
##   object - (xsubset)
##   full - (logical = FALSE)
##   size - (integer = NULL)
##   ... - ignored
##
## Rval: (character[]) variable names
##
variable.names.xsubset <- function (object, full = FALSE, size = NULL, ...) {
  ## extract all variable names
  x.names <- variable.names(object$lm)

  ## trivial
  if (full) {
    return (x.names)
  }

  ## process size
  size <- process.size(size, size.allow = object$size)

  ## compute indices
  x.which <- object$which[[which(size == object$size)]]

  ## done
  x.names[x.which]
}


## extract formula
##
## Args:
##   x - (xsubset)
##   size - (integer)
##   ... - ignored
##
## Rval: (formula)
##
formula.xsubset <- function (x, size, ...) {
  ## extract variable names
  y.name <- as.character(formula(x$lm)[[2]])
  x.names <- variable.names(x, size = size)

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
##   ... - passed to 'model.frame.lm'
##
## Rval: (data.frame)
##
model.frame.xsubset <- function (formula, ...) {
  model.frame(formula$lm, ...)
}


## extract model matrix
##
## Args:
##   object - (xsubset)
##   size - (integer)
##   ... - ignored
##
## Rval: (matrix)
##
## Weights:  The returned matrix is not affected
## by the weights (if present).
##
model.matrix.xsubset <- function (object, size, ...) {
  ## process size
  size <- process.size(size, size.allow = object$size)

  ## extract model matrix
  x <- model.matrix(object$lm)

  ## compute index
  which <- object$which[[which(size == object$size)]]

  ## done
  x[, which, drop = FALSE]
}


## extract response variable
##
## Args:
##   object - (xsubset)
##   ... - passed to 'model.response.lm'
##
## Rval: (numeric[])
##
model.response.xsubset <- function (object, ...) {
  model.response(object$lm, ...)
}


## extract model weights
##
## Args:
##   object - (xsubset)
##   ... - passed to 'model.weights.lm'
##
## Rval: (numeric[])
##
model.weights.xsubset <- function (object, ...) {
  model.weights(object$lm, ...)
}


## extract model offset
##
## Args:
##   object - (xsubset)
##   ... - passed to 'model.offset.lm'
##
## Rval: (numeric[])
##
model.offset.xsubset <- function (object, ...) {
  model.offset(object$lm, ...)
}


## extract terms object
##
## Args:
##   x - (xsubset)
##   size - (integer = NULL)
##   ... - passed to 'terms.formula'
##
## Rval: (terms)
##
## If 'size' is "NULL", then simply return the 'terms'
## of the associated 'lm' object.
##
terms.xsubset <- function (x, size = NULL, ...) {
  ## trivial
  if (is.null(size)) {
    return (terms(x$lm))
  }

  ## extract formula and model frame
  f <- formula(x, size)
  mf <- model.frame(x)

  ## done
  terms(f, data = mf, ...)
}


## refit associated 'lm' object
##
## Args:
##   object - (xsubset)
##   size - (integer)
##   ... - ignored
##   mask.call - (logical) used internally
##
## Rval: (lm)
##
## Note:  'mask.call'
## Reconstruct 'call' object for pretty printing.
##
refit.xsubset <- function (object, size, ..., mask.call = TRUE) {
  ## extract formula
  f <- formula(object, size)

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
                 data = model.frame(object),
                 weights = model.weights(object),
                 offset = model.offset(object))
    do.call("lm", args)
  }
}


## extract ceofficients
##
## Args:
##   object - (xsubset)
##   size - (integer)
##   ... - passed to 'coef.lm'
##
## Rval: (numeric[])
##
coef.xsubset <- function (object, size, ...) {
  coef(refit(object, size, mask.call = FALSE), ...)
}


## extract variance-covariance matrix
##
## Args:
##   object - (xsubset)
##   size - (integer)
##   ... - passed to 'vcov.lm'
##
## Rval: (matrix)
##
vcov.xsubset <- function (object, size, ...) {
  vcov(refit(object, size, mask.call = FALSE), ...)
}


## extract fitted values
##
## Args:
##   object - (xsubset)
##   size - (integer)
##   ... - passed to 'fitted.lm'
##
## Rval: (numeric[])
##
fitted.xsubset <- function (object, size, ...) {
  fitted(refit(object, size, mask.call = FALSE), ...)
}


## extract residuals
##
## Args:
##   object - (xsubset)
##   size - (integer)
##   ... - passed to 'residuals.lm'
##
## Rval: (numeric[])
##
residuals.xsubset <- function (object, size, ...) {
  residuals(refit(object, size, mask.call = FALSE), ...)
}


## extract deviance
##
## Args:
##   object - (xsubset)
##   size - (integer[] = NULL)
##   ... - ignored
##
## Rval: (numeric[])
##
## Returns the (weighted) RSS for the specified
## subset sizes.
##
## Note:  'size'
## The default behavior is to return the deviance
## for all selected subset sizes.
##
deviance.xsubset <- function (object, size = NULL, ...) {
  ## process size
  size <- process.size(size, size.default = object$size, size.allow = object$size,
                       allow.vector = TRUE)
  
  # compute index
  rss.which <- which(object$size %in% size)

  # done
  object$rss[rss.which]
}


## extract log-likelihood
##
## Args:
##   object - (xsubset)
##   size - (integer)
##   ... - passed to 'logLik.lm'
##
## Rval: (logLik)
##
logLik.xsubset <- function (object, size, ...) {
  ##
  ## TODO:
  ## Are the weights (if present) correctly taken
  ## into account by 'logLik.lm'?
  ##
  logLik(refit(object, size, mask.call = FALSE), ...)
}


## compute AIC
##
## Args:
##   object - (xsubset)
##   size - (integer[] = NULL)
##   k - (integer)
##   ... - ignored
##
## Rval: (numeric|data.frame)
##
## Note:  'size'
## The default behavior is to extract the AIC for
## all selected subset sizes.
##
AIC.xsubset <- function (object, size = NULL, ..., k = 2) {
  ## process size
  size <- process.size(size, size.default = object$size, size.allow = object$size,
                       allow.vector = TRUE)

  ## refits
  fitted <- lapply(size, function (sz) {
    refit(object, sz, mask.call = FALSE)
  })

  # call AIC
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
## summary
##


## summary for 'xsubset' objects
##
## Args:
##   object - (xsubset)
##   size - (integer = NULL)
##   aic.penalty - (numeric = 2)
##   ... - ignored
##
## Rval: (summary.xsubset)
##
summary.xsubset <- function (object, size = NULL, aic.penalty = 2, ...) {
  paste <- function (...) base::paste(..., sep = "")

  ## validate size
  size <- process.size(size, size.default = object$size,
                       size.allow = object$size, allow.vector = TRUE)

  ## collect information
  nreg <- object$nreg
  include <- object$include
  exclude <- object$exclude

  ## extract criteria and determine best model
  rss <- deviance(object, size)
  if (length(size) > 1) {
    aic <- AIC(object, size, k = aic.penalty)$AIC
    best.which <- which(aic == min(aic))
  } else {
    aic <- AIC(object, size, k = aic.penalty)
    best.which <- 1
  }
  best <- size[best.which]
  names(rss) <- names(aic) <- size

  ## construct summary table
  which.tab <- matrix(FALSE, ncol = length(size), nrow = nreg)
  for (i in 1:length(size)) {
    sz <- size[i]
    sz.name <- as.character(sz)
    sz.which <- object$which[[sz.name]]
    which.tab[sz.which, i] <- TRUE
  }
  rownames(which.tab) <- variable.names(object, full = TRUE)
  colnames(which.tab) <- size
  colnames(which.tab)[best.which] <- paste(colnames(which.tab)[best.which], "*")
  rownames(which.tab)[include] <- paste("+", rownames(which.tab)[include])
  rownames(which.tab)[exclude] <- paste("-", rownames(which.tab)[exclude])
  
  ## replace some slots with updated information
  object$size <- size
  object$which <- which.tab
  object$rss <- rss
  object$aic <- aic
  object$best <- best
  object$lm <- summary(object$lm)
  
  ## object class
  class(object) <- "summary.xsubset"

  ## done
  object
}


## print 'xsubset' summary
##
## Args:
##   x - (summary.xsubset)
##   digits - (integer = ...)
##   ... - ignored
##
## Rval: (summary.xsubset) invisible
##
print.summary.xsubset <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
  catln <- function (...) base::cat(..., "\n", sep = "")

  catln()
  catln("Call:")
  catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))

  catln()
  catln("Selected variables:")
  which.x <- ifelse(x$which, "x", "")
  print(which.x, quote = FALSE)

  catln()
  catln("Model fit:")
  fit <- rbind(x$aic, x$rss)
  rownames(fit) <- c("AIC", "RSS")
  colnames(fit) <- colnames(which.x)
  print(fit, digits = digits)
    
  invisible(x)
}


## plot 'xsubset' summary
##
## Args:
##   x - (xsubset)
##   size - (integer = NULL)
##   legend - (logical = TRUE)
##   xlab - (character = ...)
##   ylab - (character = "")
##   main - (character = NULL)
##   col - (integer[]|character[] = ...)
##   lty - (integer = 1)
##   type - (character = "b")
##   ... -
##
## Rval: (summary.xsubset) invisible
##
## All arguments are passed to 'plot.default'.
##
plot.summary.xsubset <- function (x, legend = TRUE,
  xlab = "Number of regressors in model", ylab = "",
  main = NULL, col = c("blue", "red"), lty = 1, type = "b", ...)
{
  ## extract main info
  size <- x$size
  rss <- x$rss
  aic <- x$aic

  col <- rep(col, length.out = 2)
  type <- rep(type, length.out = 2)
  lty <- rep(lty, length.out = 2)
  if (is.null(main)) {
    main <- paste("AIC and residual sum of squares")
  }

  plot(size, rss, ylab = ylab, xlab = xlab, main = main,
    type = type[1], lty = lty[1], col = col[1], ...)
  new.old <- getOption("new")
  par(new = TRUE)
  plot(size, aic, type = type[2], axes = FALSE, col = col[2],
       xlab = "", ylab = "", ...)
  if (legend) {
    legend("topright", c("RSS", "AIC"),
           lty = lty, col = col, bty = "n")
  }
  axis(4)
  par(new = new.old)
  
  invisible(x)
}


##
## Misc
##


process.size <- function (size, size.default = NULL, size.allow = NULL,
                          allow.vector = FALSE, canonicalize = FALSE) {
  if (missing(size)) {
    stop ("missing argument 'size'")
  }
  if (is.null(size)) {
    size <- size.default
  }
  if (!is.numeric(size)) {
    stop ("'size' must be numeric")
  }
  if (!allow.vector && (length(size) > 1)) {
    size <- size[1]
    warning ("redundant values in 'size'; using ", size)
  }
  if (!is.null(size.allow) && !all(size %in% size.allow)) {
    stop ("requested 'size' not available")
  }
  size.rounded <- round(size)
  if (any(size != size.rounded)) {
    warning ("'size' must be integer; rounding 'size'")
    size <- size.rounded
  }
  if (canonicalize) {
    size <- sort(unique(size))
  }
  size
}


pretty.integer <- function (val) {
  if (length(val) < 1) {
    return ("-")
  }
  val.min <- min(val)
  val.max <- max(val)
  if (val.min == val.max) {
    format(val)
  } else if (all(val == val.min:val.max)) {
    paste(val.min, val.max, sep = ":")
  } else if (is.vector(val)) {
    paste(val, collapse = ", ")
  } else {
    "-"
  }
}
