## new generics
xsubset <- function(object, ...) UseMethod("xsubset")

refit <- function(object, ...) UseMethod("refit")


## xsubset methods
## standard formula interface
xsubset.formula <- function(formula, data, subset, na.action, weights, offset, ...,
  model = TRUE, y = TRUE, x = FALSE)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## call model.frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract X and Y
  mt <- terms(formula, data = data)
  X <- model.matrix(mt, mf)
  Y <- model.extract(mf, "response")

  ## weights and offset
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- rep(1, NROW(X))
  weights <- structure(as.vector(weights), .Names = rownames(mf))
  offset <- model.offset(mf)
  if(is.null(offset)) offset <- rep(0, NROW(X))
  offset <- as.vector(offset)

  ## call default method
  rval <- xsubset.default(X, Y, weights, offset, na.action = na.action, ...)

  ## add formula/model.frame type information
  rval$call <- call
  rval$formula <- formula
  rval$terms <- mt
  rval$model <- if(model) mf else NULL
  if(!y) rval$y <- NULL
  if(!x) rval$x <- NULL
  if(!is.null(attr(mf, "na.action"))) rval$na.action <- attr(mf, "na.action")

  return(rval)
}

## interface for fitted lm regression
xsubset.lm <- function(object, ..., model = TRUE, y = TRUE, x = FALSE)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## extract information from fitted model object
  mt <- terms(object)
  mf <- model.frame(object)
  X <- if(is.matrix(object$x)) object$x else model.matrix(mt, mf)
  Y <- if(is.vector(object$y)) object$y else model.response(mf)
  weights <- weights(object)
  offset <- object$offset

  ## call default method
  rval <- xsubset.default(X, Y, weights, offset, ...)
  
  ## add formula/model.frame type information
  rval$call <- call
  rval$formula <- formula(object)
  rval$terms <- mt
  rval$model <- if(model) mf else NULL
  if(!y) rval$y <- NULL
  if(!x) rval$x <- NULL
  if(!is.null(attr(mf, "na.action"))) rval$na.action <- attr(mf, "na.action")  
  rval$contrasts <- object$contrasts
  
  return(rval)
}

## workhorse method
xsubset.default <- function(x, y, weights = NULL, offset = NULL,
  size = NULL, pradius = 1, tol = 0,
  na.action = na.omit)
#Z# Argument | Previously | Comment
#Z# ---------+------------+--------------------------------------------------------
#Z# size     | vmin/vmax  | -
#Z# pradius  | prad       | Shouldn't the default be NROW(x)/3?
#Z# tol      | tau        | What is the range of this? Could we have tol in [0, 1]?
#Z# weights  | -          | Could weights (for WLS rather than OLS) be incorporated
#Z#          |            | into the underlying C code?
#Z# offset   | -          | Additional offset: E[y] = X beta + offset
#Z#          |            | Can easily be handled by internal transformation
#Z#          |            | of y (done).
#Z# ?        | -          | Is there an argument that can force a particular
#Z#          |            | variable/column in or out of the regression?
#Z#          |            | Maybe 'mark' does something like that?
#Z# ?        | -          | Is there some possibility to group columns together,
#Z#          |            | i.e., include all or none? (could be used for factors)
{
  ## keep call (of generic)
  call <- match.call()
  call[[1]] <- as.name("xsubset")

  ## process x, y, weights, offset
  x <- as.matrix(x)
  stopifnot(is.numeric(x), is.numeric(y), NROW(x) == NROW(y))
  if(is.null(offset)) offset <- rep(0, length(y))  
  stopifnot(length(offset) == length(y))
  ay <- cbind(x, y - offset)
  if(any(is.na(ay))) {
    ay <- na.action(ay)
    na.action <- attr(ay, "na.action")
  } else {
    na.action <- NULL
  }
  nobs <- nrow(ay)
  nvar <- ncol(ay)
  if(is.null(weights)) weights <- rep(1, nobs)
  if(!isTRUE(all.equal(as.vector(weights), rep(1, nobs)))) stop("'weights' not implemented yet")
  #Z# Could weights be added to the underlying C code?

  ## intercept processing
  icpt <- isTRUE(all.equal(as.vector(ay[,1]), rep(1, nobs)))
  if(icpt) {
    size <- if(is.null(size)) 1:(nvar-1) else sort(unique(round(size))) + 1
    mark <- 1 #Z# What exactly is 'mark'?
    null.rss <- var(y - offset) * (nobs - 1)
  } else {
    size <- if(is.null(size)) 0:(nvar-1) else sort(unique(round(size)))
    mark <- 0
    null.rss <- sum((y - offset)^2)
  }
  if(max(size) > nvar - 1) stop("'size' can not be larger than number of regressors")
  if(min(size) < icpt) stop("'size' must at least be 0")
  if(as.integer(icpt) %in% size) {
    null <- TRUE
    size <- size[-which(size == as.integer(icpt))]
  } else {
    null <- FALSE
  }
  
  ## tolerance
  .TAU_MAX <- .Machine$double.xmax #Z# Can we have tol in [0, 1] or some other fixed interval?
  ntol <- rep(.TAU_MAX, length.out = nvar - 1)
  ntol[size] <- rep(tol, length.out = length(size))
  tol <- ntol
  
  ## call underlying C code
  C_rval <- .C(name = "xsubset_R",
    ## in
    nobs = as.integer(nobs),
    nvar = as.integer(nvar),
    ay   = as.numeric(ay),
    mark = as.integer(mark),
    prad = as.integer(pradius),
    tau  = as.numeric(tol),
    ## out
    rsel = numeric(nvar - 1),
    isel = integer(nvar * (nvar - 1) / 2),
    nvis = integer(1)
  )

  ## extract selected variables and associated RSS
  vwhich <- sapply(size, function(i) {
    j <- i * (i - 1) / 2
    C_rval$isel[(j + 1 + icpt):(j + i)] + 1
  })
  rss <- C_rval$rsel[size]

  ## include null?
  if(null) {
    size <- c(as.integer(icpt), size)
    vwhich <- c(list("0" = NULL), vwhich)
    rss <- c(null.rss, rss)
  }
  
  ## switch back to size without intercept
  size <- size - icpt
  names(vwhich) <- names(rss) <- size

  ## compute BIC and best model
  logL <- -0.5 * nobs * (log(rss) + 1 - log(nobs) + log(2 * pi))
  bic <- -2 * logL + log(nobs) * (size + icpt + 1)
  names(logL) <- names(bic) <- size
  best <- size[which.min(bic)]

  rval <- list(
    call = call,
    size = size,
    rss = rss,
    which = vwhich,
    bic = bic,
    best = best,
    x = x,
    y = y,
    weights = NULL, #Z# if supported: if(isTRUE(all.equal(as.vector(weights), rep(1, nobs)))) NULL else weights
    offset = if(isTRUE(all.equal(as.vector(offset), rep(0, nobs)))) NULL else offset,
    nobs = nobs,
    nreg = nvar - 1 - icpt,
    intercept = icpt,
    pradius = pradius,
    tol = tol,
    nvis = C_rval$nvis, #Z# What is 'nvis'?
    na.action = na.action
  )
  class(rval) <- "xsubset"
  return(rval)
}

print.xsubset <- function(x, ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), sep = "\n")

  ## format size
  size <- x$size
  size <- if(isTRUE(all.equal(as.vector(size), min(size):max(size)))) paste(min(size), max(size), sep = ":")
    else paste(size, collapse = ", ")
  
  rval <- format(c(x$nobs, x$nreg, x$best), width = max(nchar(size), 3))
  rval <- as.matrix(c(rval[1:2],
    format(if(x$intercept) "yes" else "no", width = max(nchar(size)), justify = "right"),
    size, rval[3]))
  colnames(rval) <- ""
  rownames(rval) <- c("Number of observations:", "Number of regressors:",
    "Intercept:", "Subset sizes assessed:", "Best BIC subset size:")
  print(rval, quote = FALSE)
  #Z# Should we also include technical arguments pradius/tol?
  
  invisible(x)
}

summary.xsubset <- function(object, size = NULL, ...) {

  ## get size and associated which
  if(is.null(size)) size <- object$size
  which <- object$which  

  ## get maximal model matrix
  x <- if(is.null(object$x)) model.matrix(terms(object), model.frame(object)) else object$x

  ## collect information
  rss <- deviance(object, size = size)
  bic <- AIC(object, size = size, k = log(object$nobs))
  wbest <- which(size == object$best)
  nwhich <- matrix(FALSE, ncol = length(size), nrow = ncol(x))
  for(i in 1:length(size)) {
    wi <- which[[as.character(size[i])]]
    if(object$intercept) wi <- c(1, wi)
    if(!is.null(wi)) nwhich[wi, i] <- TRUE
  }
  rownames(nwhich) <- colnames(x)
  colnames(nwhich) <- size
  colnames(nwhich)[wbest] <- paste(colnames(nwhich)[wbest], "*", sep = "")
  
  ## replace some slots with updated information
  object$size <- size
  object$which <- nwhich
  object$rss <- rss
  object$bic <- bic

  ## delete some slots
  object$x <- object$y <- object$formula <- object$terms <- object$model <-
    object$na.action <- object$contrasts <- object$weights <- object$offset <- NULL
  class(object) <- "summary.xsubset"

  return(object)
}

print.summary.xsubset <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  cat("Selected variables:\n")
  wi <- x$which
  wi[] <- ifelse(wi, "x", "")
  print(wi, quote = FALSE)

  cat("\nModel fit:\n")
  fit <- rbind(x$bic, x$rss)
  rownames(fit) <- c("BIC", "RSS")
  colnames(fit) <- colnames(wi)
  print(fit, digits = digits)
    
  invisible(x)
}

plot.xsubset <- function(x, legend = TRUE,
  xlab = "Number of regressors in model", ylab = "",
  main = "BIC and residual sum of squares",
  col = c(1, 4), lty = 1, type = "b", ...)
{
  ## extract main info
  size <- x$size
  rss <- deviance(x, size = size)
  bic <- AIC(x, size = size, k = log(x$nobs))

  col <- rep(col, length.out = 2)
  type <- rep(type, length.out = 2)
  lty <- rep(lty, length.out = 2)

  plot(size, bic,
    ylab = ylab, xlab = xlab, main = main,
    type = type[1], lty = lty[1], col = col[1], ...)
  onew <- getOption("new")
  par(new = TRUE)
  plot(size, rss, type = type[2], axes = FALSE, col = col[2], xlab = "", ylab = "")
  if(legend) legend("topright", c("BIC", "RSS"), lty = lty, col = col, bty = "n")
  axis(4)
  par(new = onew)
  invisible(x)
}

deviance.xsubset <- function(object, size = NULL, ...) {
  rsize <- if(is.null(size)) object$best else size
  size <- object$size
  if(!all(rsize %in% size)) stop("requested 'size' not available")
  rval <- object$rss[which(size %in% rsize)]
  if(length(rval) < 2) rval <- as.vector(rval)
  return(rval)
}

logLik.xsubset <- function(object, size = NULL, ...) {
  size <- if(is.null(size)) object$best else size[1]
  rss <- deviance(object, size = size)
  structure(-0.5 * object$nobs * (log(rss) + 1 - log(object$nobs) + log(2 * pi)),
    df = size + object$intercept + 1, class = "logLik")
}

AIC.xsubset <- function(object, size = NULL, ..., k = 2) {
  rss <- deviance(object, size = size)
  logl <- -0.5 * object$nobs * (log(rss) + 1 - log(object$nobs) + log(2 * pi))  
  structure(-2 * logl + k * (size + object$intercept + 1), .Names = names(rss))
}

model.frame.xsubset <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  if(is.null(formula$formula)) stop("'model.frame' not available")
  NextMethod()
}

model.matrix.xsubset <- function(object, size = NULL, ...) {
  size <- if(is.null(size)) object$best else size[1]
  if(!any(size == object$size)) stop("requested 'size' not available")
  x <- object$x
  if(is.null(object$x)) x <- model.matrix(terms(object), model.frame(object), contrasts = object$contrasts)
  wi <- which(size == object$size)
  wi <- object$which[[wi]]
  if(object$intercept) wi <- c(1, wi)
  if(is.null(wi)) wi <- 0
  x[, wi, drop = FALSE]
}

refit.xsubset <- function(object, size = NULL, ...) {
  size <- if(is.null(size)) object$best else size[1]

  ## extract x/y
  x <- model.matrix(object, size = size)
  if(object$intercept) x <- x[, -1, drop = FALSE]
  y <- if(!is.null(object$y)) object$y else model.response(model.frame(object))
  
  ## set up data.frame and formula
  data <- cbind(data.frame(y = y), as.data.frame(x))
  if(!is.null(object$formula)) names(data)[1] <- as.character(object$formula[[2]])  
  formula <- as.formula(paste(names(data)[1], "~",
    if(!object$intercept) "0 +" else NULL, 
    paste(names(data)[-1], collapse = " + ")))
  
  ## offset and weights
  offset <- object$offset
  weights <- object$weights
  
  lm(formula, data, weights = weights, offset = offset)
}

coef.xsubset <- function(object, size = NULL, ...) {
  coef(refit(object, size), ...)
}

vcov.xsubset <- function(object, size = NULL, ...) {
  vcov(refit(object, size), ...)
}

