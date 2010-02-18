##
## File:  xselect.R
##
## Author:  Marc Hofmann
##
## Based on:  xsubset.R
##
## :-|  so much for code reuse...
##


## method for matrix objects
##
## Args:
##   object  - (matrix) model matrix
##   y       - (~matrix = NULL) response variable
##   ... - arguments forwarded to 'xselect.data.frame'
##
## Rval: (list)
##   (see 'xselect.default')
##
## This method coerces the matrix 'cbind(y, x)' to a
## 'data.fame' and forwards the call to 'xselect.data.frame'.
##
## Note:  response variable
## If 'y' is 'NULL', the last column of 'x' is used as the
## response variable.
##
## Note:  intercept
## No automatic intercept detection.
## 
##
xselect.matrix <- function (object, y = NULL, ...) {
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xselect")

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
  rval <- xselect(df, ...)

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
##   ... - arguments forwarded to 'xselect.formula'
##
## Rval: (list)
##   (see 'xselect.default')
##
## This method simply constructs a 'formula' from the
## object ('formula(object, env = NULL)') and forwards
## the call to 'xselect.formula'.
##
xselect.data.frame <- function (object, y = 1, ...) {
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xselect")

  ## rearrange data frame
  if (is.character(y)) {
    y <- match(y, names(object))
  }
  df <- data.frame(object[y], object[-y])

  ## extract formula
  f <- formula(df, env = NULL)

  ## forward call
  rval <- xselect(f, object, ...)

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
##   ... - arguments forwarded to 'xselect.lm'
##
## Rval: (list)
##   (see 'xselect.default')
##
## Constructs an 'lm' object and forwards the call
## to 'xselect.lm'.
##
## Note:  arguments
## The named arguments are passed to 'lm'.
##
xselect.formula <- function (object, data = NULL, row.subset = NULL, weights = NULL,
                             na.action = na.omit, model = TRUE, x = FALSE, y = FALSE,
                             contrasts = NULL, offset = NULL, ...)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xselect")

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
  rval <- xselect(lm, ...)

  ## update return value
  rval$call <- call

  ## done
  rval
}


## interface for fitted lm regression
##
## Args:
##   lm - (lm)
##   ... - arguments forwarded to 'xselect.default'
##
## Rval:
##   (see 'xselect.default')
##
## Does nothing.  Simply forwards to 'xselect.default'.
##   
xselect.lm <- function (object, ...)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xselect")

  ## forward call
  rval <- xselect.default(object, ...)

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
##   tolerance - (numeric[] = 0)
##   criterion - (numeric|character)
##   ... - ignored
##   pradius - (integer)
##   precision - (character)
##
## Rval: (list)
##   call - (call)
##   lm   - (lm)
##   nreg - (integer)
##   include - (integer[])
##   exclude - (integer[])
##   criterion - (character|numeric)
##   which - (list)
##   value - (numeric)
##   intercept - (logical)
##
## All the fun happens here.
##
xselect.default <- function (object, include = NULL, exclude = NULL,
                             tolerance = 0, criterion = NULL, ...,
                             pradius = NULL, precision = "double")
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xselect")

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

  ## process tolerance
  tolerance.computed <- tolerance
  if (!is.numeric(tolerance.computed)) {
    stop ("'tolerance' must be numeric")
  }
  if (length(tolerance.computed) > 1) {
    tolerance.computed <- tolerance.computed[1]
    warning ("redundant values in 'tolerance'; using ", tolerance.computed)
  }

  ## process criterion
  if (is.null(criterion)) {
    criterion <- "AIC"
  }
  criterion.computed <- criterion
  penalty.computed <- 0
  if (is.numeric(criterion.computed)) {
    criterion.computed <- "mAIC"
    penalty.computed <- criterion
  } else if (!is.character(criterion.computed)) {
    stop ("'criterion' must be numeric or character")    
  }
  criterion.computed <- tolower(criterion.computed)
  if (!(criterion.computed %in% c("aic", "aicc", "aicc2", "aicu",
                                  "bic", "cp", "maic"))) {
    stop ("unknown 'criterion'")
  }

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
                 criterion = as.character(criterion.computed),
                 penalty   = as.double(penalty.computed),
                 nobs      = as.integer(nobs.computed),
                 nreg      = as.integer(nreg.computed),
                 xy        = as.double(cbind(x.computed, y.computed)),
                 mark      = as.integer(mark.computed),
                 tolerance = as.double(tolerance.computed),
                 pradius   = as.integer(pradius.computed),
                 ## out
                 value = double(1),
                 which = as.integer(rep(-1, nreg.computed)),
                 count = integer(1))
  C_rval <- do.call(".C", c(name = "R_select", C_args))

  ## extract selected regressors
  which <- C_rval$which + 1
  which <- column.index[which[which > 0]]
  best <- length(which)

  ## extract criterion value
  value <- C_rval$value
  
  ## names
  names(which) <- names(value) <- best

  # return value
  rval <- list(call      = call,
               lm        = object,
               nreg      = nreg,
               include   = include.computed,
               exclude   = exclude.computed,
               criterion = criterion,
               best      = best,
               which     = which,
               value     = value,
               intercept = intercept.computed,
               .nodes = C_rval$count)
  class(rval) <- "xselect"

  # done
  rval
}


## print method for 'xselect' objects
##
## Args:
##   x - (xselect)
##   ... - ignored
##
## Rval: (xselect) invisible
##
print.xselect <- function (x, ...)
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

  criterion.pretty <- pretty.criterion(x$criterion)
  best.pretty <- format(x$best)

  val <- as.matrix(c(criterion.pretty, best.pretty), ncol = 1)
  colnames(val) <- ""
  rownames(val) <- c("Criterion:", "Best subset:")
  print(val, quote = FALSE)

  invisible(x)
}


## extract environment
##
## Args:
##   object - (xselect)
##   ... - passed to 'environment.default'
##
## Rval: (environment)
##
## Returns the 'environment' object associated with
## this object's 'lm' object.
##
environment.xselect <- function (object, ...) {
  environment(object$lm, ...)
}


## extract variable names
##
## Args:
##   object - (xselect)
##   ... - passed to 'variable.names.lm'
##
## Rval: (character[]) variable names
##
variable.names.xselect <- function (object, full = FALSE, ...) {
  ## extract all variable names
  x.names <- variable.names(object$lm, ...)

  ## trivial
  if (full) {
    ## done
    return (x.names)
  }

  ## done
  x.names[object$which]
}


## extract formula
##
## Args:
##   x - (xselect)
##   ... - ignored'
##
## Rval: (formula)
##
formula.xselect <- function (x, ...) {
  ## extract variable names
  y.name <- as.character(formula(x$lm)[[2]])
  x.names <- variable.names(x)

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
##   formula - (xselect)
##   ... - passed to 'model.frame.lm'
##
## Rval: (data.frame)
##
model.frame.xselect <- function (formula, ...) {
  model.frame(formula$lm, ...)
}


## extract model matrix
##
## Args:
##   object - (xselect)
##   ... - passed to 'model.matrix.lm'
##
## Rval: (matrix)
##
model.matrix.xselect <- function (object, ...) {
  ## extract model matrix
  x <- model.matrix(object$lm, ...)

  ## done
  x[, object$which, drop = FALSE]
}


## extract response variable
##
## Args:
##   object - (xselect)
##   ... - passed to 'model.response.lm'
##
## Rval: (numeric[])
##
model.response.xselect <- function (object, ...) {
  model.response(object$lm, ...)
}


## extract model weights
##
## Args:
##   object - (xselect)
##   ... - passed to 'model.weights.lm'
##
## Rval: (numeric[])
##
model.weights.xselect <- function (object, ...) {
  model.weights(object$lm, ...)
}


## extract model offset
##
## Args:
##   object - (xselect)
##   ... - passed to 'model.offset.lm'
##
## Rval: (numeric[])
##
model.offset.xselect <- function (object, ...) {
  model.offset(object$lm, ...)
}


## extract terms object
##
## Args:
##   x - (xselect)
##   ... - passed to 'terms.formula'
##
## Rval: (terms)
##
terms.xselect <- function (x, ...) {
  ## extract formula and model frame
  f <- formula(x)
  mf <- model.frame(x)

  ## done
  terms(f, data = mf, ...)
}


## refit associated 'lm' object
##
## Args:
##   object - (xselect)
##   ... - ignored
##   mask.call - (logical)
##
## Rval: (lm)
##
refit.xselect <- function (object, ..., mask.call = TRUE) {
  ## extract formula
  f <- formula(object)

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
##   object - (xselect)
##   ... - passed to 'coef.lm'
##
## Rval: (numeric[])
##
coef.xselect <- function (object, ...) {
  coef(refit(object, mask.call = FALSE), ...)
}


## extract variance-covariance matrix
##
## Args:
##   object - (xselect)
##   ... - passed to 'vcov.lm'
##
## Rval: (matrix)
##
vcov.xselect <- function (object, ...) {
  vcov(refit(object, mask.call = FALSE), ...)
}


## extract fitted values
##
## Args:
##   object - (xselect)
##   ... - passed to 'fitted.lm'
##
## Rval: (numeric[])
##
fitted.xselect <- function (object, ...) {
  fitted(refit(object, mask.call = FALSE), ...)
}


## extract residuals
##
## Args:
##   object - (xselect)
##   ... - passed to 'residuals.lm'
##
## Rval: (numeric)
##
residuals.xselect <- function (object, ...) {
  residuals(refit(object, mask.call = FALSE), ...)
}


## extract deviance
##
## Args:
##   object - (xselect)
##   ... - passed to 'deviance.lm'
##
## Rval: (numeric)
##
deviance.xselect <- function (object, ...) {
  deviance(refit(object, mask.call = FALSE), ...)
}


## extract log-likelihood
##
## Args:
##   object - (xselect)
##   ... - passed to 'logLik.lm'
##
## Rval: (logLik)
##
logLik.xselect <- function (object, ...) {
  logLik(refit(object, mask.call = FALSE), ...)
}


## compute AIC
##
## Args:
##   object - (xselect)
##   k - (integer)
##   ... - passed to 'AIC.default'
##
## Rval: (numeric)
##
AIC.xselect <- function (object, ..., k = 2) {
  AIC(refit(object, mask.call = FALSE, k = k), ...)
}


##
## summary
##


## summary for 'xselect' objects
##
## Args:
##   object - (xselect)
##   ... - ignored
##
## Rval: (summary.xselect)
##
summary.xselect <- function (object, ...) {
  paste <- function (...) base::paste(..., sep = "")
  
  ## collect information
  nreg <- object$nreg
  include <- object$include
  exclude <- object$exclude
  best <- object$best
  which <- object$which

  ## collect information
  which.tab <- (1:nreg) %in% which
  which.tab <- matrix(which.tab, ncol = 1, nrow = nreg)
  rownames(which.tab) <- variable.names(object, full = TRUE)
  colnames(which.tab) <- best
  rownames(which.tab)[include] <- paste("+", rownames(which.tab)[include])
  rownames(which.tab)[exclude] <- paste("-", rownames(which.tab)[exclude])
  
  ## update information
  object$which <- which.tab
  object$lm <- summary(object$lm)

  ## class
  class(object) <- "summary.xselect"

  object
}


## print 'xselect' summary
##
## Args:
##   x - (summary.xselect)
##   digits - (integer = ...)
##   ... - ignored
##
## Rval: (summary.xselect) invisible
##
print.summary.xselect <- function (x, digits = max(3, getOption("digits") - 3), ...)
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
  fit <- rbind(x$value)
  rownames(fit) <- pretty.criterion(x$criterion)
  colnames(fit) <- x$best
  print(fit, digits = digits)
    
  invisible(x)
}


##
## Misc
##

pretty.criterion <- function (x, ...) {
  paste <- function (...) base::paste(..., sep = "")

  if (is.numeric(x)) {
    return (paste("mAIC (penalty = ", format(x), ")"))
  }
  x <- tolower(x)
  if (x == "aic") {
    "AIC"
  } else if (x == "aicc") {
    "AICc"
  } else if (x == "aicc2") {
    "AICc2"
  } else if (x == "aicu") {
    "AICu"
  } else if (x == "bic") {
    "BIC"
  } else if (x == "cp") {
    "Cp"
  } else {
    stop ("undefined 'criterion'")
  }
}
