##
## File:  lmSubsets.R
##



###########################
##  WORKHORSE FUNCTIONS  ##
###########################


## workhorse function
##
## Args:
##   x         - (matrix) model matrix
##   y         - (numeric[]) model response
##   weights   - (numeric[]) 
##   offset    - (numeric[])
##   include   - (integer[]) regressors to force in
##   exclude   - (integer[]) regressors to force out
##   size      - (integer[]) subset sizes
##   tolerance - (numeric[])
##   pradius   - (integer) preordering radius
##   nbest     - (integer) number of best subsets
##   ...       - ignored
##   .algo     - (character) for interal use
##
## Rval:  (list)
##   nobs      - (integer)
##   nvar      - (integer)
##   weights   - (numeric[]) 
##   offset    - (numeric[])
##   intercept - (logical)
##   include   - (integer[])
##   exclude   - (integer[])
##   size      - (integer[])
##   tolerance - (numeric[])
##   nbest     - (integer)
##   df        - (integer[nbest,nvar])
##   rss       - (numeric[nbest,nvar])
##   which     - (numeric[nvar,nbest,nvar])
##   .info     - (integer)
##   .nodes    - (integer)
##
## NOTE: '.algo'
##   The following values are recognized:
##   'dca'  :  Dropping Column Algorithm (no branch and bound)
##   'bba'  :  Branch and Bound Algorithm
##   'pbba' :  BBA with variable preordering
##   'hbba' :  Heuristic BBA
##   'hpbba':  Heuristic, preordering BBA
##   'xbba.':  experimental PBBA
##
lmSubsets.fit <- function (x, y, weights = NULL, offset = NULL,
                           include = NULL, exclude = NULL, size = NULL,
                           tolerance = 0, pradius = NULL, nbest = 1, ...,
                           .algo = "hpbba") {
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
    x.names <- sapply(colnames(x), as.name)

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
            warning ("invalid columns selected in 'exclude'; fixing 'exclude': ",
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
    has.intercept <- all(x[, 1] == 1)
    if (has.intercept) {
        if (any(include == 1)) {
            ## OK, already selected
        } else if (any(exclude == 1)) {
            ## not selected
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
        pmin <- size.max - size.min - pradius + 1
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
    C_df <- rep(-1, nvar * nbest)
    C_rss <- rep(0, nvar * nbest)
    C_which <- rep(0, nvar * nbest * nvar)

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
        df      = as.integer(C_df),
        rss     = as.double(C_rss),
        which   = as.logical(C_which),
        tau     = as.double(C_tau),
        info    = integer(1),
        nodes   = integer(1),
        NAOK    = TRUE
    )

    C_rval <- do.call(".C", c("R_lmSubsets", C_args))
    if (C_rval$info < 0) {
        stop ("invalid argument: ", names(C_rval)[-C_rval$info])
    }

    ## return value
    rval <- list(call      = call,
                 nobs      = nobs,
                 nvar      = nvar,
                 weights   = weights,
                 intercept = has.intercept,
                 include   = include,
                 exclude   = exclude,
                 size      = size,
                 tolerance = NULL,
                 nbest     = nbest,
                 df        = NULL,
                 rss       = NULL,
                 which     = NULL,
                 .info     = C_rval$info,
                 .nodes    = C_rval$nodes)
    class(rval) <- "lmSubsets"

    ## extract value and subsets
    .cnames <- format.ordinal(1:nbest)
    ## tolerance
    rval$tolerance <- array(tolerance, dim = nvar,
                            dimnames = list(1:nvar))
    ## df
    rval$df <- array(C_rval$df, dim = c(nbest, nvar),
                     dimnames = list(.cnames, 1:nvar))
    not.valid <- rval$df < 0
    rval$df[not.valid] <- NA
    ## rss
    rval$rss <- array(C_rval$rss, dim = c(nbest, nvar),
                      dimnames = list(.cnames, 1:nvar))
    rval$rss[not.valid] <- NA
    ## which
    rval$which <- array(C_rval$which, dim = c(nvar, nbest, nvar),
                        dimnames = list(x.names, .cnames, 1:nvar))
    rval$which[rep(not.valid, each = nvar)] <- NA

    ## done
    rval
}



##########################
##  GENERATOR FUNCTION  ##
##########################


## standard formula interface
##
## Args:
##   formula   - (formula) an object that can be coerced to class 'formula'
##   data      - (data.frame)
##   subset    - (vector) subset of observations
##   weights   - (numeric[])
##   na.action - (function)
##   model     - (logical)
##   x         - (logical)
##   y         - (logical)
##   contrasts - (list)
##   offset    - (numeric[])
##   ...       - forwarded to 'lmSubsets.fit'
##
## Rval:  (lmSubsets) see 'lmSubsets.fit'
##
lmSubsets <- function (formula, data, subset, weights, na.action,
                       model = TRUE, x = FALSE, y = FALSE,
                       contrasts = NULL, offset, ...) {
    ret.m <- model;  model <- NULL
    ret.x <- x;  x <- NULL
    ret.y <- y;  y <- NULL

    ## keep call (of generic)
    cl <- match.call()
    cl[[1]] <- as.name("lmSubsets")

    ## model frame
    mf <- match.call(expand.dots = FALSE)
    m <- c("formula", "data", "subset",
           "weights", "na.action", "offset")
    m <- match(m, names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    ## model terms
    mt <- attr(mf, "terms")

    ## model matrix and response
    x <- model.matrix(mt, mf, contrasts)
    y <- model.response(mf, "numeric")

    ## weights
    w <- as.vector(model.weights(mf))
    if(!is.null(w) && !is.numeric(w)) {
        stop ("'weights' must be a numeric vector")
    }

    ## offset
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset) && (length(offset) != NROW(y))) {
        stop ("length of 'offset' must equal number of observations")
    }

    ## fit subsets
    rval <- lmSubsets.fit(x, y, weights = w, offset = offset, ...)

    ## return value
    class(rval) <- "lmSubsets"
    rval$call <- cl
    rval$na.action <- attr(mf, "na.action")
    rval$offset <- offset
    rval$contrasts <- attr(x, "contrasts")
    rval$xlevels <- .getXlevels(mt, mf)
    rval$terms <- mt

    if (ret.m) rval$model <- mf
    if (ret.x) rval$x <- x
    if (ret.y) rval$y <- y

    ## done
    rval
}



########################
##  STANDARD METHODS  ##
########################


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

    x.names <- dimnames(x$which)[[1]]


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
                           c("N observations", "N regressors", "Weights", "Offset",
                             "Intercept", "Include", "Exclude", "Subset sizes",
                             "Tolerance", "N best"),
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


    ## which
    catln()
    catln("Which (size):")
    catln("  variable x best")
    which <- apply(x$which[, , , drop = FALSE], c(1, 2), format.which)
    rownames(which) <- paste("  ", rownames(which))
    print(which, quote = FALSE)


    catln()

    ## done
    invisible(x)
}


## plot method for 'lmSubsets' objects
##
## Args:
##   x    -   (lmSubsets)
##   ...  -   passed to plotting functions
##   legend - (character)
##
## Rval:  (lmSubsets) invisible
##
plot.lmSubsets <- function (x, ..., legend) {
    localPlot <- function (object, main, sub = NULL, xlab, ylab, type = "o",
                           lty = c(1, 3), pch = c(16, 21), col = "black",
                           bg = "white", ...) {
        if (missing(main)) main <- "All subsets"
        if (missing(xlab)) xlab <- "Number of regressors"
        if (missing(ylab)) ylab <- "Deviance"

        type <- rep(type, length = 2)
        lty  <- rep(lty, length = 2)
        pch  <- rep(pch, length = 2)
        col  <- rep(col, length = 2)
        bg   <- rep(bg, length = 2)

        x <- matrix(rep(object$size, each = object$nbest), nrow = object$nbest)
        y <- object$rss[object$nbest:1, object$size, drop = FALSE]

        matplot(x = x, y = y, main = main, sub = sub, xlab = xlab, ylab = ylab,
                type = type[2], lty = lty[2], pch = pch[2], col = col[2],
                bg = bg[2], ...)
        lines(x[1, ], y[object$nbest, ], type = type[1], lty = lty[1],
              pch = pch[1], col = col[1], bg = bg[1], ...)

        if (!is.null(legend)) {
            legend("topright", legend = legend, lty = lty,
                   pch = pch, col = col, pt.bg = bg, bty = "n")
        }
    }

    if (inherits(x, "summary.lmSubsets")) {
        if (missing(legend)) legend <- "Deviance (RSS)"

        localPlot(x, ...)
    } else {
        plot(summary(x, ...))
    }

    invisible(x)
}



#########################
##  EXTRACTOR METHODS  ##
#########################


## extract variable names
##
## Args:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (character[])
##   variable names
##
variable.names.lmSubsets <- function (object, size, best = 1, ...) {
    x.names <- dimnames(object$which)[[1]]

    if (missing(size)) {
        return (x.names)
    }

    x.names[object$which[, best, size]]
}


## extract formula
##
## Args:
##   x    - (lmSubsets)
##   size - (integer)
##   best - (integer)
##   ...  - ignored
##
## Rval:  (formula)
##
formula.lmSubsets <- function (x, size, best = 1, ...) {
    if (missing(size)) {
        f <- formula(x$terms)

        return (f)
    }

    e <- new.env();
    e$x <- model.matrix(x, size = size, best = best)
    e$y <- model.response(x)

    if (x$intercept) {
        e$x <- e$x[, -1]
        f <- formula("y ~ x + 1", env = e)
    } else {
        f <- formula("y ~ x + 0", env = e)
    }
    
    f
}


## extract model frame
##
## Args:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (character[])
##   variable names
##
model.frame.lmSubsets <- function(formula, size, best = 1, ...) {
    mf <- formula[["model"]]

    if (!missing(size) || is.null(mf)) {
        cl <- formula$call
        m <- c("formula", "data", "subset",
               "weights", "na.action", "offset")
        m <- match(m, names(cl), 0L)
        cl <- cl[c(1L, m)]
        cl$drop.unused.levels <- TRUE
        cl[[1L]] <- quote(stats::model.frame)
        cl$xlev <- formula$xlevels
        cl$formula <- stats::formula(formula, size = size, best = best)

        env <- environment(formula$terms)
        if (is.null(env)) {
            env <- parent.frame()
        }

        mf <- eval(cl, env)
    }

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
    x <- object[["x"]]
    if (is.null(x)) {
        data <- model.frame(object)
        x <- model.matrix.default(object, data = data, contrasts.arg = object$contrasts)
    }

    if (missing(size)) {
        return (x)
    }

    x[, object$which[, best, size]]
}


## extract model response
##
## Args:
##   data - (lmSubsets)
##   ...  - ignored
##
## Rval:  (formula)
##
model.response.lmSubsets <- function (data, ...) {
    y <- data[["y"]]
    if (is.null(y)) {
        mf <- model.frame(data)
        y <- model.response(mf)
    }

    y
}


## refit
##
## Args:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - forwarded to 'lm'
##
## Rval:  (lm)
##
refit.lmSubsets <- function (object, size, best = 1, ...) {
    f <- formula(object, size = size, best = best)

    lm(f, ...)
}


## extract ceofficients
##
## Args:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   This method refits the 'lmSubsets' object and calls 'coef' on the
##   obtained 'lm' object.  In some circumstances it might be more
##   efficient to call 'refit' directly and execute 'coef' on the
##   obtained 'lm' object.
##
coef.lmSubsets <- function (object, size, best = 1, ...) {
    coef(refit(object, size = size, best = best, model = FALSE))
}


## extract variance-covariance matrix
##
## Args:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (matrix)
##
## Note: 'refit'
##   See 'coef.lmSubsets'.
##
vcov.lmSubsets <- function (object, size, best = 1, ...) {
    vcov(refit(object, size = size, best = best, model = FALSE))
}


## extract fitted values
##
## Args:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.lmSubsets'.
##
fitted.lmSubsets <- function (object, size, best = 1, ...) {
    fitted(refit(object, size = size, best = best, model = FALSE))
}


## extract residuals
##
## Args:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.lmSubsets'.
##
residuals.lmSubsets <- function (object, size, best = 1, ...) {
    residuals(refit(object, size = size, best = best, model = FALSE))
}


## extract deviance
##
## Args:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (numeric[])
##
## Returns the RSS for the specified subset.
##
deviance.lmSubsets <- function (object, size, best = 1, ...) {
    ## size
    if (missing(size)) {
        d <- deviance(refit(object, model = FALSE))

        return (d)
    }

    ## extract RSS
    object$rss[best, size]
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
logLik.lmSubsets <- function (object, size, best = 1, ...) {
    ## size
    if (missing(size)) {
        ll <- logLik(refit(object, model = FALSE))

        return (ll)
    }

    ## weights
    N <- object$nobs
    if (is.null(w <- object$weights)) {
        S <- 0
    } else {
        S <- sum(log(w[w != 0]))
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
