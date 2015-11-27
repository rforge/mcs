##
## File:  lmSubsets.R
##



################
## GENERATORS ##
################


## standard formula interface
##
## Args:
##   formula - (formula)
##   ...     - forwarded to 'lm' and 'lmSubsets.lm'
##   lm      - (logical) if 'true', compute 'lm' component
##
## Rval:  (lmSubsets)
##   See 'lmSubsets.default'.
##
## NOTE:  'lm'
##   If 'lm' is 'FALSE', the returned 'lm' component is
##   an "empty" mockup.
##
lmSubsets.formula <- function (formula, ..., lm = FALSE) {
    ret.lm <- lm;  lm <- NULL

    ## keep call (of generic)
    call <- match.call()
    call[[1]] <- as.name("lmSubsets")

    ## lm call
    lm.call <- match.call(expand.dots = TRUE)
    m <- match(c("include", "exclude", "size", "tolerance",
                 "pradius", "nbest", "lm"), names(lm.call), 0L)
    lm.call[[1]] <- as.name("lm")
    lm.call[m] <- NULL

    ## build 'lm' component
    if (ret.lm) {
        object <- eval(lm.call, parent.frame())

        ret.m <- TRUE
        ret.x <- TRUE
        ret.y <- TRUE
    } else {
        ## mock-up
        dots <- list(...)
        if (is.null(ret.m <- dots[["model"]])) ret.m <- TRUE
        if (is.null(ret.x <- dots[["x"]])) ret.x <- FALSE
        if (is.null(ret.y <- dots[["y"]])) ret.y <- FALSE
        ## model frame
        mf.call <- match.call(expand.dots = TRUE)
        m <- match(c("formula", "data", "subset",
                     "weights", "na.action", "offset"),
                   names(mf.call), 0L)
        mf.call <- mf.call[c(1L, m)]
        mf.call$drop.unused.levels <- TRUE
        mf.call[[1L]] <- quote(stats::model.frame)
        mf <- eval(mf.call, parent.frame())
        ## model terms
        mt <- attr(mf, "terms")
        ## model matrix and response
        x <- model.matrix(mt, mf, dots$contrasts)
        y <- model.response(mf, "numeric")
        ## mock fit
        qr <- array(rep(NA, NCOL(x)), dim = c(1, NCOL(x)),
                    dimnames = list(NULL, colnames(x)))
        ## mock lm
        object <- list(call = lm.call)
        class(object) <- "lm"
        object$rank <- NCOL(x)
        object$na.action <- attr(mf, "na.action")
        object$offset <- as.vector(model.offset(mf))
        object$contrasts <- attr(x, "contrasts")
        object$xlevels <- .getXlevels(mt, mf)
        object$terms <- mt
        object$model <- mf
        object$x <- x
        object$y <- y
        object$qr <- list(qr = qr)
    }

    ## forward call
    cl <- match.call(expand.dots = TRUE)
    cl[[1L]] <- as.name("lmSubsets")
    cl$formula <- cl$y <- cl$lm <- NULL
    cl$object <- object
    rval <- eval(cl, parent.frame())

    ## return value
    rval$call <- call
    if (!ret.m) rval$lm$model <- NULL
    if (!ret.x) rval$lm$x <- NULL
    if (!ret.y) rval$lm$y <- NULL
    if (!ret.lm) rval$lm <- NULL

    ## done
    rval
}


## interface for fitted lm regressions
##
## Args:
##   object - (lm)
##   ...    - forwarded to 'lmSubsets.default'
##
## Rval:  (lmSubsets)
##   See 'lmSubsets.default'.
##
lmSubsets.lm <- function (object, ...) {
    ## keep call (of generic)
    call <- match.call()
    call[[1]] <- as.name("lmSubsets")

    ## model frame and terms
    mf <- model.frame(object)
    mt <- terms(object)
    if (is.empty.model(mt)) {
        stop ("empty model")
    }
    if (any(attr(mt, "dataClasses") == "factor")) {
        stop ("model contains factors")
    }

    ## model matrix and response
    x <- model.matrix(object)
    y <- model.response(mf)

    ## model weights and offset
    w <- model.weights(mf)
    o <- model.offset(mf)

    ## forward call
    rval <- lmSubsets(x, y, weights = w, offset = o, ...)

    ## return value
    rval$call <- call
    rval$.lm <- rval$lm <- object

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
##   tolerance - 
##   pradius   - (integer) preordering radius
##   nbest     - (integer) number of best subsets
##   ...       - ignored
##   .algo     - (character) for interal use
##
## Rval:  (lmSubsets)
##   call      - (call)
##   nobs      - (integer)
##   nvar      - (integer)
##   weights   - (numeric[])
##   offset    - (numeric[])
##   include   - (integer[])
##   exclude   - (integer[])
##   size      - (integer[])
##   intercept - (logical)
##   tolerance - (numeric[])
##   nbest     - (integer)
##   df        - (integer[nbest,nvar])
##   rss       - (numeric[nbest,nvar])
##   which     - (numeric[nvar,nbest,nvar])
##   .nodes    - (integer)
##
## NOTE: '.algo'
##   The following values are recognized:
##   'dca'  :  Dropping Column Algorithm (no branch and bound)
##   'bba'  :  Branch and Bound Algorithm
##   'pbba' :  BBA with variable preordering
##   'hbba' :  PBBA with tolerances
##   'xbba*':  experimental PBBA
##
lmSubsets.default <- function (object, y, weights = NULL, offset = NULL,
                               include = NULL, exclude = NULL,
                               size = NULL, tolerance = 0, pradius = NULL,
                               nbest = 1, ..., .algo = "hbba") {
    ## keep call (of generic)
    call <- match.call()
    call[[1]] <- as.name("lmSubsets")

    ## model matrix
    x <- as.matrix(object);  object <- NULL

    ## intercept
    intercept <- all(x[, 1] == 1)

    ## model weights and offset
    if (is.null(w <- weights)) w <- rep(1, NROW(x))
    if (is.null(o <- offset)) o <- rep(0, NROW(x))
    wok <- w != 0;  w <- w[wok];  o <- o[wok]
    x <- sqrt(w) * x[wok, , drop = FALSE]
    y <- sqrt(w) * y[wok,   drop = FALSE] - o

    ## dims
    nobs <- NROW(x)
    nvar <- NCOL(x)

    ## variables names
    x.names <- colnames(x)

    ## include
    if (is.null(include)) {
        include <- numeric(0)
    } else {
        ## numeric --> numeric
        if (is.numeric(include)) {
            include <- match(include, 1:nvar)
        }
        ## character --> numeric
        if (is.character(include)) {
            include <- match(include, x.names)
        }
        ## logical --> numeric
        if (is.logical(include)) {
            include <- which(rep(include, length.out = nvar))
        }
        ## canonicalize include
        if (any(is.na(include))) {
            include <- na.omit(include)
            warning ("non-existing columns selected in 'include'; fixing 'include': ",
                     "c(", paste(include, collapse = ","), ")")
        }
        include <- sort(unique(include))
    }

    ## exclude
    if (is.null(exclude)) {
        exclude <- numeric(0)
    } else {
        ## numeric --> numeric
        if (is.numeric(exclude)) {
            exclude <- match(exclude, 1:nvar)
        }
        ## character --> numeric
        if (is.character(exclude)) {
            exclude <- match(exclude, x.names)
        }
        ## logical --> numeric
        if (is.logical(exclude)) {
            exclude <- which(rep(exclude, length.out = nvar))
        }
        ## canonicalize exclude
        if (any(is.na(include))) {
            exclude <- na.omit(exclude)
            warning ("invalid columns selected in 'excluded'; fixing 'exclude': ",
                     "c(", paste(exclude, collapse = ","), ")")
        }
        exclude <- sort(unique(exclude))
    }

    ## include/exclude non-overlapping
    if (any(intersect(include, exclude))) {
        exclude <- setdiff(exclude, include)
        warning ("'include' and 'exclude' overlap; fixing 'exclude': ",
                 "c(", paste(exclude, collapse = ","), ")")
    }

    ## intercept
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

    ## which
    which <- c(include, setdiff(1:nvar, include))
    which <- setdiff(which, exclude)

    ## size
    size.min <- length(include) + 1
    size.max <- length(which)
    if (is.null(size)) {
        size <- size.min:size.max
    } else {
        size <- match(size, size.min:size.max)
        size <- (size.min:size.max)[size]
        if (any(is.na(size))) {
            size <- na.omit(size)
            warning ("invalid sizes selected in 'size'; fixing 'size': ",
                     "c(", paste(size, collapse = ","), ")")
        }
    }

    ## preordering radius
    ## TODO: validate
    if (is.null(pradius)) {
        pmin <- 8
    } else {
        pmin <- size.max - pradius
    }

    ## tolerance
    ## TODO: validate
    tolerance <- rep(tolerance, length.out = nvar)
    tolerance[-size] <- +Inf

    ## call C
    C_mark <- size.min - 1
    C_size <- size.max
    C_v <- which - 1
    C_xy <- cbind(x[, which], y)
    C_tau <- tolerance + 1

    C_index <- rep(0, nvar * nbest)
    C_rss <- rep(0, nvar * nbest)
    C_which <- rep(0, nvar * nbest * nvar)

    C_tau <- C_tau[size]
    C_tau[C_tau == +Inf] <- .Machine$double.xmax
    C_tau <- c(0, C_tau)

    C_args <- list(
        algo    = as.character(.algo),
        nobs    = as.integer(nobs),
        nvar    = as.integer(nvar),
        size    = as.integer(C_size),
        mark    = as.integer(C_mark),
        nbest   = as.integer(nbest),
        pmin    = as.integer(pmin),
        v       = as.integer(C_v),
        xy      = as.double(C_xy),
        rss     = as.double(C_rss),
        which   = as.logical(C_which),
        tau     = as.double(C_tau),
        nodes   = integer(1),
        info    = integer(1)
    )

    C_rval <- do.call(".C", c("R_lmSubsets", C_args))

    ## return value
    rval <- list(call      = call,
                 nobs      = nobs,
                 nvar      = nvar,
                 weights   = weights[wok],
                 offset    = offset[wok],
                 include   = include,
                 exclude   = exclude,
                 size      = NULL,
                 intercept = intercept,
                 tolerance = NULL,
                 nbest     = nbest,
                 df        = NULL,
                 rss       = NULL,
                 which     = NULL,
                 .nodes = C_rval$nodes)
    class(rval) <- "lmSubsets"

    ## extract value & subsets
    .cnames <- paste(1:nbest, ".", sep = "")
    ## size
    rval$size <- size
    ## tolerance
    rval$tolerance <- array(tolerance, dim = nvar,
                            dimnames = list(1:nvar))
    ## df
    rval$df <- array(rep(1:nvar, each = nbest) + 1, dim = c(nbest, nvar),
                     dimnames = list(.cnames, 1:nvar))
    ## rss
    rval$rss <- array(C_rval$rss, dim = c(nbest, nvar),
                      dimnames = list(.cnames, 1:nvar))
    rval$rss[, -rval$size] <- NA
    rval$rss[rval$rss == .Machine$double.xmax] <- NA
    ## which
    rval$which <- array(C_rval$which, dim = c(nvar, nbest, nvar),
                        dimnames = list(x.names, .cnames, 1:nvar))
    rval$which[, , -rval$size] <- NA

    ## done
    rval
}


## print method for 'lmSubsets' objects
##
## Args:
##   x   - (lmSubsets)
##   ... - ignored
##
## Rval:  (lmSubsets) invisible
##
print.lmSubsets <- function (x, ...) {
    catln <- function (...) base::cat(..., "\n", sep = "")
    paste <- function (..., sep = "") base::paste(..., sep = sep)

    x.names <- variable.names(x, .full = TRUE)


    ## call
    catln()
    catln("Call:")
    catln("  ", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))


    ## arguments
    catln()
    cat("Arguments:")
    val <- as.matrix(c(format(x$nobs),
                       format(x$nvar),
                       if (is.null(x$weights)) "No" else "Yes",
                       if (is.null(x$offset)) "No" else "Yes",
                       if (x$intercept) "Yes" else "No",
                       paste(x.names[x$include], collapse = " "),
                       paste(x.names[x$exclude], collapse = " "),
                       paste(format(rbind(x$size, format(x$tolerance[x$size], nsmall = 2)))[1, ], collapse = "|"),
                       paste(format(rbind(x$size, format(x$tolerance[x$size], nsmall = 2)))[2, ], collapse = "|"),
                       format(x$nbest)),
                     ncol = 1)
    colnames(val) <- ""
    rownames(val) <- paste("  ",
                           c("Number of observations", "Number of regressors",
                             "Weights", "Offset", "Intercept", "Include", "Exclude",
                             "Subset sizes", "Tolerance", "N best"),
                           ":")
    print(val, quote = FALSE)

    ## fit
    catln()
    catln("Model fit (deviance):")
    catln("  best x size")
    rss <- x$rss[, x$size, drop = FALSE]
    fit <- ifelse(is.na(rss), "", format(rss, nsmall = 2))
    rownames(fit) <- paste("  ", rownames(fit))
    print(fit, quote = FALSE)
    catln()

    ## done
    invisible(x)
}


## plot method for 'lmSubsets' objects
##
## Args:
##   x   - (lmSubsets)
##   ... - ignored
##   add - (logical) add to current plot
##
## Rval:  (lmSubsets) invisible
##
## Graphical arguments are passed to 'plot.default'.
##
plot.lmSubsets <- function (x, ..., add = FALSE) {
    if (!add) {
        xlim <- range(x$size)
        ylim <- range(x$rss[!is.na(x$rss)])

        ## start
        plot.new()
        box()

        ## title
        title(main = "lmSubsets")

        ## legend
        legend("topright", legend = "Deviance",
               lty = 3, pch = 21, pt.bg = "white")
        
        ## window
        plot.window(xlim = xlim, ylim = ylim)
        axis(1, at = pretty(xlim))
        axis(2, at = pretty(ylim))
        mtext("Size", side = 1, line = 3, col = "black")
        mtext("Deviance", side = 2, line = 3, col = "black")
    }

    ## plot deviance
    matplot(x = matrix(rep(x$size, each = x$nbest), nrow = x$nbest),
            y = x$rss[, x$size, drop = FALSE],
            type = "o", lty = 3, pch = 21,
            col = "black", bg = "white",
            add = TRUE)

    ## done
    invisible(x)
}



#############
## METHODS ##
#############


## extract variable names
##
## Args:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##   .full  - (logical) for internal use
##   .cmpl  - (logical) for internal use
##
## Rval:  (character[])
##   variable names
##
variable.names.lmSubsets <- function (object, size, best = 1, ...,
                                      .full = FALSE, .cmpl = FALSE) {
    ## lm
    if (is.null(object$.lm)) {
        stop ("'lmSubsets' object does not have a 'lm' component")
    }

    ## full model
    x.names <- variable.names(object$.lm)
    if (.full) {
        return (x.names)
    }

    ## which
    wi <- object$which[, best, size]
    ## complementary model
    wi <- xor(wi, .cmpl)
    ## variable names
    x.names <- x.names[wi]

    ## done
    x.names
}


## extract formula
##
## Args:
##   x   - (lmSubsets)
##   ... - forwarded to 'variable.names.lmSubsets'
##
## Rval:  (formula)
##
formula.lmSubsets <- function (x, ...) {
    ## lm
    if (is.null(x$.lm)) {
        stop ("'lmSubsets' object does not have a 'lm' component")
    }

    ## full model
    f <- formula(x$.lm)
    ## variable names (complementary submodel)
    x.names <- variable.names(x, ..., .cmpl = TRUE)
    ## update formula
    f <- update(f, paste(".~.", paste(c("", x.names), collapse = "-")))

    ## done
    f
}


## extract model frame
##
## Args:
##   formula - (lmSubsets)
##   ...     - forwarded to 'formula.lmSubsets'
##
## Rval:  (data.frame)
##
model.frame.lmSubsets <- function (formula, ...) {
    ## lm
    if (is.null(formula$.lm)) {
        stop ("'lmSubsets' object does not have a 'lm' component")
    }

    ## full model
    mf <- model.frame(formula$.lm)
    ## formula
    f <- formula(formula, ...)
    ## build
    mf <- model.frame(f, data = mf)

    ## done
    mf
}


## extract model matrix
##
## Args:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (matrix)
##
model.matrix.lmSubsets <- function (object, size, best = 1, ...) {
    ## lm
    if (is.null(object$.lm)) {
        stop ("'lmSubsets' object does not have a 'lm' component")
    }

    ## full model
    x <- model.matrix(object$.lm)
    ## which
    wi <- object$which[, best, size]
    ## select
    x <- x[, wi]

    ## done
    x
}


## refit
##
## Args:
##   object - (lmSubsets)
##   ...    - forwarded to 'formula.lmSubsets'
##
## Rval:  (lm)
##
refit.lmSubsets <- function (object, ...) {
    ## lm
    if (is.null(object$.lm)) {
        stop ("'lmSubsets' object does not have a 'lm' component")
    }

    ## extract formula
    f <- formula(object, ...)
    ## model frame call
    mf.call <- match.call(expand.dots = TRUE)
    mf.call[[1]] <- as.name("model.frame")
    mf.call$formula <- mf.call$object;  mf.call$object <- NULL
    ## lm call
    lm.call <- object$.lm$call
    lm.call$formula <- f
    lm.call$data <- mf.call
    ## eval
    lm <- eval(lm.call, parent.frame())

    ## done
    lm
}


## extract ceofficients
##
## Args:
##   object - (lmSubsets)
##   ...    - forwarded to 'refit.lmSubsets'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   This method refits the 'lmSubsets' object and calls 'coef' on the
##   obtained 'lm' object.  In some circumstances it might be more
##   efficient to call 'refit' directly and execute 'coef' on the
##   obtained 'lm' object.
##
coef.lmSubsets <- function (object, ...) {
    coef(refit(object, ...))
}


## extract variance-covariance matrix
##
## Args:
##   object - (lmSubsets)
##   ...    - forwarded to 'refit.lmSubsets'
##
## Rval:  (matrix)
##
## Note: 'refit'
##   See 'coef.lmSubsets'.
##
vcov.lmSubsets <- function (object, ...) {
    vcov(refit(object, ...))
}


## extract fitted values
##
## Args:
##   object - (lmSubsets)
##   ...    - forwarded to 'refit.lmSubsets'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.lmSubsets'.
##
fitted.lmSubsets <- function (object, ...) {
  fitted(refit(object, ...))
}


## extract residuals
##
## Args:
##   object - (lmSubsets)
##   ...    - forwarded to 'refit.lmSubsets'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.lmSubsets'.
##
residuals.lmSubsets <- function (object, ...) {
  residuals(refit(object, ...))
}


## extract deviance
##
## Args:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (numeric[])
##
## Returns the RSS for the specified subset(s).
##
deviance.lmSubsets <- function (object, size, best = 1, ...) {
    ## size
    if (missing(size)) size <- object$size

    ## extract RSS
    d <- object$rss[best, size]

    ## done
    d
}


## extract log-likelihood (base method)
##
## Args:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer)
##   ...    - ignored
##   df     - (integer[]) degrees of freedom
##
## Rval: (logLik)
##
## Note:
##   Can handle multiple objects.
##
logLik.lmSubsets <- function (object, size, best = 1, ...) {
    ## size
    if (missing(size)) size <- object$size

    ## weights
    N <- object$nobs
    if (is.null(w <- object$weights)) {
        S <- 0
    } else {
        S <- sum(log(w))
    }

    ## extract rss
    rss <- deviance(object, size = size, best = best)

    ## done
    structure(0.5 * (S - N * (log(2 * pi) + 1 - log(N) + log(rss))),
              df       = object$df[best, size]                     ,
              nobs     = N                                         ,
              class    = "logLik"                                  )
}


## compute AIC
##
## Args:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer)
##   ...    - ignored
##   k      - (numeric)
##
## Rval: (numeric|data.frame)
##
AIC.lmSubsets <- function (object, size, best = 1, ..., k = 2) {
    ## size
    if (missing(size)) size <- object$size

    ## extract log-likelihoods
    ll <- logLik(object, size = size, best = best)
    ## compute AICs
    aic <- AIC(ll, k = k)
    ## data frame?
    if (length(aic) > 1) {
        aic <- data.frame(df  = attr(ll, "df"),
                          AIC = aic           )
    }

    ## done
    aic
}


## compute BIC
##
## Args:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer)
##   ...    - ignored
##
## Rval: (numeric|data.frame)
##
BIC.lmSubsets <- function (object, size, best = 1, ...) {
    ## size
    if (missing(size)) size <- object$size

    ## extract log-likelihoods
    ll <- logLik(object, size = size, best = best)
    ## compute BICs
    bic <- BIC(ll)
    ## data frame?
    if (length(bic) > 1) {
        bic <- data.frame(df  = attr(ll, "df"),
                          BIC = bic           )
    }

    ## done
    bic
}
