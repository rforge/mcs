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
##   nmin      - (integer)
##   nmax      - (integer)
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
##   nmin      - (integer)
##   nmax      - (integer)
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
##   'phbba':  HBBA with variable preordering
##   'xbba.':  experimental PBBA
##
lmSubsets_fit <- function (x, y, weights = NULL, offset = NULL,
                           include = NULL, exclude = NULL, nmin = NULL,
                           nmax = NULL, tolerance = 0, pradius = NULL,
                           nbest = 1, ..., .algo = "phbba") {
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
    size <- length(which)
    mark <- length(include)
    nmin <- max(nmin, mark + 1)
    nmax <- min(nmax, size)

    ## preordering radius
    ## TODO: validate
    if (is.null(pradius)) {
        pmin <- 8
    } else {
        pmin <- size - mark - pradius + 1
    }

    ## tolerance
    ## TODO: validate
    tolerance <- rep(tolerance, length.out = nvar)
    tolerance[-(nmin:nmax)] <- +Inf

    ## call C
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
        size    = as.integer(size),
        mark    = as.integer(mark),
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
                 nmin      = nmin,
                 nmax      = nmax,
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



#########################
##  GENERATOR METHODS  ##
#########################


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
##   ...       - forwarded to 'lmSubsets_fit'
##
## Rval:  (lmSubsets) see 'lmSubsets_fit'
##
lmSubsets.default <- function (formula, data, subset, weights, na.action,
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
    rval <- lmSubsets_fit(x, y, weights = w, offset = offset, ...)

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

    ## call
    catln()
    catln("Call:")
    catln("  ", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))

    ## arguments
    catln()
    cat("Arguments:")
    output <- c(format(x$nobs), format(x$nvar),
                if (is.null(x$weights)) "No" else "Yes",
                if (is.null(x$offset)) "No" else "Yes",
                if (x$intercept) "Yes" else "No",
                format(x$nmin), format(x$nmax),
                format(x$nbest))
    output <- as.matrix(output, ncol = 1)
    colnames(output) <- ""
    rownames(output) <- paste("  ", c("N observations", "N regressors",
                                      "Weights", "Offset", "Intercept",
                                      "N min", "N max", "N best"), ":")
    print(output, quote = FALSE)

    ## fit
    catln()
    catln("Model fit:")
    catln("  best x size (tolerance)")
    catln("    DEVIANCE")
    output <- x$rss[, x$nmin:x$nmax, drop = FALSE]
    output <- ifelse(is.na(output), "", format(output, nsmall = 2))
    rownames(output) <- paste("  ", rownames(output))
    colnames(output) <- paste(colnames(output), " (", x$tolerance[x$nmin:x$nmax], ")")
    print(output, quote = FALSE)

    ## which
    catln()
    catln("Which:")
    catln("  variable x best")
    catln("    Size")
    output <- apply(x$which[, , , drop = FALSE], c(1, 2), format.which)
    rownames(output)[x$include] <- paste("+", rownames(output)[x$include])
    rownames(output)[x$exclude] <- paste("-", rownames(output)[x$exclude])
    rownames(output) <- paste("  ", rownames(output))
    print(output, quote = FALSE)

    ## pad
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
    localPlot <- function (object, main, sub, xlab, ylab, type = "o",
                           lty = c(1, 3), pch = c(16, 21), col = "black",
                           bg = "white", ...) {
        if (missing(main)) main <- "All subsets"
        if (missing(sub))  sub  <- paste("nbest =", object$nbest, sep = " ")
        if (missing(xlab)) xlab <- "Number of regressors"
        if (missing(ylab)) ylab <- "Deviance"

        type <- rep(type, length = 2)
        lty  <- rep(lty, length = 2)
        pch  <- rep(pch, length = 2)
        col  <- rep(col, length = 2)
        bg   <- rep(bg, length = 2)

        x <- matrix(rep(object$nmin:object$nmax, each = object$nbest),
                    nrow = object$nbest)
        y <- object$rss[object$nbest:1, object$nmin:object$nmax, drop = FALSE]

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

    ## default legend
    if (missing(legend)) legend <- "Deviance (RSS)"

    ## plot
    localPlot(x, ...)

    ## done
    invisible(x)
}



#########################
##  EXTRACTOR METHODS  ##
#########################


## extract variable names
##
## Args:
##   object - (lmSubsets)
##   size   - (integer|character)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (character[])
##
variable.names.lmSubsets <- function (object, size, best = 1, ...) {
    ## full model
    x.names <- dimnames(object$which)[[1]]

    ## 'size' processing
    if (missing(size)) {
        return (x.names)
    }

    ## shortcut to 'lmSelect'
    if (is.character(size)) {
        if (!missing(best))  warning ("unused argument: 'best'")

        size <- lmSelect(object, penalty = size)$df[1] - 1
    }

    ## submodel
    which <- object$which[, best, size]
    x.names[which]
}


## extract formula
##
## Args:
##   x    - (lmSubsets)
##   size - (integer|character)
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

    ## shortcut to 'lmSelect'
    if (is.character(size)) {
        if (!missing(best))  warning ("unused argument: 'best'")

        size <- lmSelect(object, penalty = size)$df[1] - 1
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
##   size   - (integer|character)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (model.frame)
##
model.frame.lmSubsets <- function (formula, size, best = 1, ...) {
    mf <- formula[["model"]]

    ## shortcut to 'lmSelect'
    if (!missing(size) && is.character(size)) {
        if (!missing(best))  warning ("unused argument: 'best'")

        size <- lmSelect(object, penalty = size)$df[1] - 1
    }

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
##   object        - (lmSubsets)
##   size          - (integer|character)
##   best          - (integer)
##   ...           - ignored
##
## Rval:  (numeric[,])
##
model.matrix.lmSubsets <- function (object, size, best = 1, ...) {
    ## full model
    x <- object[["x"]]
    if (is.null(x)) {
        data <- model.frame(object)
        x <- model.matrix.default(object, data = data,
                                  contrasts.arg = object$contrasts)
    }

    ## 'size' processing
    if (missing(size)) {
        return (x)
    }

    ## shortcut to 'lmSelect'
    if (is.character(size)) {
        if (!missing(best))  warning ("unused argument: 'best'")

        size <- lmSelect(object, penalty = size)$df[1] - 1
    }

    ## submodel
    which <- object$which[, best, size]
    x[, which]
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

    ## done
    y
}


## refit
##
## Args:
##   object - (lmSubsets)
##   size   - (integer|character)
##   best   - (integer)
##   ...    - forwarded to 'lm'
##
## Rval:  (lm)
##
refit.lmSubsets <- function (object, size, best = 1, ...) {
    if (missing (size)) {
        stop ("missing argument: 'size'")
    }

    ## shortcut to 'lmSelect'
    if (is.character(size)) {
        if (!missing(best))  warning ("unused argument: 'best'")

        size <- lmSelect(object, penalty = size)$df[1] - 1
    }

    f <- formula(object, size = size, best = best)

    lm(f, ...)
}


## extract deviance
##
## Args:
##   object - (lmSubsets)
##   size   - (integer[]|character)
##   best   - (integer[])
##   ...    - ignored
##   drop   - (logical)
##
## Rval:  (numeric[,])
##
deviance.lmSubsets <- function (object, size, best = 1, ..., drop = TRUE) {
    ## 'size' processing
    if (missing(size)) {
        size <- object$nmin:object$nmax
    } else if (is.character(size)) {
        if (!missing(best))  warning ("unused argument: 'best'")

        size <- lmSelect(object, penalty = size)$df[1] - 1
    }

    ## extract RSS
    ans <- object$rss[best, size, drop = FALSE]
    ans <- t(ans)

    ## drop
    if (drop) {
        ans <- drop(ans)
    }

    ## done
    ans
}


## extract log-likelihood (base method)
##
## Args:
##   object - (lmSubsets)
##   size   - (integer[]|character)
##   best   - (integer[])
##   ...    - ignored
##   drop   - (logical)
##
## Rval: (numeric[,])
##
logLik.lmSubsets <- function (object, size, best = 1, ..., drop = TRUE) {
    ## 'size' processing
    if (missing(size)) {
        size <- object$nmin:object$nmax
    } else if (is.character(size)) {
        if (!missing(best))  warning ("unused argument: 'best'")

        size <- lmSelect(object, penalty = size)$df[1] - 1
    }

    ## weights
    N <- object$nobs
    if (is.null(w <- object$weights)) {
        S <- 0
    } else {
        S <- sum(log(w[w != 0]))
    }

    ## extract rss
    rss <- deviance(object, size = size, best = best, drop = FALSE)
    if (drop) rss <- as.numeric(rss)

    ## extract df
    df <- object$df[best, size, drop = FALSE]
    df <- t(df)
    if (drop) df <- as.numeric(df)

    ## compute log-likelihoods
    ll <- 0.5 * (S - N * (log(2 * pi) + 1 - log(N) + log(rss)))
    ans <- structure(ll, df = df, nobs = N)
    if (drop) class(ans) <- "logLik"

    ## done
    ans
}


## compute AIC
##
## Args:
##   object - (lmSubsets)
##   size   - (integer[]|character)
##   best   - (integer[])
##   ...    - ignored
##   k      - (numeric) penalty
##   drop   - (logical)
##
## Rval: (numeric[,])
##
AIC.lmSubsets <- function (object, size, best = 1, ..., k = 2, drop = TRUE) {
    ## extract log-likelihoods
    ll <- match.call(expand.dots = FALSE)
    m <- c("object", "size", "best", "drop")
    m <- match(m, names(ll), 0L)
    ll <- ll[c(1L, m)]
    ll[[1L]] <- quote(stats::logLik)
    ll <- eval(ll)

    ## compute AICs
    ans <- -2 * as.numeric(ll) + k * attr(ll, "df")

    ## data frame?
    if (drop && (length(ans) > 1)) {
        ans <- data.frame(df = attr(ll, "df"), AIC = ans)
    }

    ## done
    ans
}


## compute BIC
##
## Args:
##   object - (lmSubsets)
##   size   - (integer[]|character)
##   best   - (integer[])
##   ...    - ignored
##   drop   - (logical)
##
## Rval: (numeric[,])
##
BIC.lmSubsets <- function (object, size, best = 1, ..., drop = TRUE) {
    ## extract log-likelihoods
    ll <- match.call(expand.dots = FALSE)
    m <- c("object", "size", "best", "drop")
    m <- match(m, names(ll), 0L)
    ll <- ll[c(1L, m)]
    ll[[1L]] <- quote(stats::logLik)
    ll <- eval(ll)

    ## compute BICs
    k <- log(attr(ll, "nobs"))
    ans <- -2 * as.numeric(ll) + k * attr(ll, "df")

    ## data frame?
    if (drop && (length(ans) > 1)) {
        ans <- data.frame(df = attr(ll, "df"), BIC = ans)
    }

    ## done
    ans
}


## extract ceofficients
##
## Args:
##   object        - (lmSubsets)
##   size          - (integer|character)
##   best          - (integer)
##   ...           - ignored
##
## Rval:  (numeric[])
##
coef.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    ## shortcut to 'lmSelect'
    if (is.character(size)) {
        if (!missing(best))  warning ("unused argument: 'best'")

        size <- lmSelect(object, penalty = size)$df[1] - 1
    }

    ## extract submodel
    x <- model.matrix(object, size = size, best = best)
    y <- model.response(object)

    ## solve
    ans <- qr.solve(x, y)

    ## names
    x.names <- variable.names(object, size = size, best = best)
    names(ans) <- x.names

    ## done
    ans
}


## extract variance-covariance matrix
##
## Args:
##   object        - (lmSubsets)
##   size          - (integer|character)
##   best          - (integer)
##   ...           - ignored
##
## Rval:  (numeric[,])
##
vcov.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    ## shortcut to 'lmSelect'
    if (is.character(size)) {
        if (!missing(best))  warning ("unused argument: 'best'")

        size <- lmSelect(object, penalty = size)$df[1] - 1
    }

    ## extract submodel
    x <- model.matrix(object, size = size, best = best)
    y <- model.response(object)

    ## compute 'vcov'
    qr <- qr(cbind(x, y))
    r <- qr$qr[1:size, 1:size, drop = FALSE]
    rdf <- object$nobs - size
    rss <- object$rss[best, size]
    sigma2 <- rss/rdf
    ans <- chol2inv(r) * sigma2

    ## names
    x.names <- variable.names(object, size = size, best = best)
    dimnames(ans) <- list(x.names, x.names)

    ## done
    ans
}


## extract fitted values
##
## Args:
##   object - (lmSubsets)
##   size   - (integer|character)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (numeric[])
##
fitted.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    ## shortcut to 'lmSelect'
    if (is.character(size)) {
        if (!missing(best))  warning ("unused argument: 'best'")

        size <- lmSelect(object, penalty = size)$df[1] - 1
    }

    ## extract submodel
    x <- model.matrix(object, size = size, best = best)
    y <- model.response(object)

    ## fit
    qr <- qr(x)
    qr.fitted(qr, y)
}


## extract residuals
##
## Args:
##   object - (lmSubsets)
##   size   - (integer|character)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (numeric[])
##
residuals.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    ## shortcut to 'lmSelect'
    if (is.character(size)) {
        if (!missing(best))  warning ("unused argument: 'best'")

        size <- lmSelect(object, penalty = size)$df[1] - 1
    }

    ## extract submodel
    x <- model.matrix(object, size = size, best = best)
    y <- model.response(object)

    ## fit
    qr.x <- qr(x)
    qr.resid(qr.x, y)
}
