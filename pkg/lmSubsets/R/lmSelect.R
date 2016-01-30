##
## File:  lmSelect.R
##



###########################
##  WORKHORSE FUNCTIONS  ##
###########################


## coerce 'lmSubsets' to 'lmSelect' object
##
## Args:
##   object   - (lmSubsets)
##   penalty  - (numeric) penalty per parameter
##   ...      - ignored
##
## Rval:  (lmSelect)
##   See 'lmSelect.default'.
##
lmSubsets.select <- function (object, penalty = "BIC", ...) {
    ## tolerance
    tolerance <- object$tolerance[object$nmin:object$nmax]
    if (!all(tolerance[1] == tolerance)) {
        stop ("non-homogeneous tolerance")
    }
    object$tolerance <- tolerance[1]

    ## call
    object$call[1] <- call("lmSelect")

    ## aic
    N <- object$nobs
    if (is.null(w <- object$weights)) {
        S <- 0
    } else {
        S <- sum(log(w[w != 0]))
    }

    if (tolower(penalty) == "aic") {
        object$penalty <- structure("AIC", k = 2)
    } else if (tolower(penalty) == "bic") {
        object$penalty <- structure("BIC", k = log(N))
    } else if (is.numeric(penalty)) {
        object$penalty <- structure("AIC", k = penalty)
    } else {
        stop ("invalid 'penalty'")
    }

    ll <- 0.5 * (S - N * (log(2 * pi) + 1 - log(N) + log(object$rss)))
    aic <- -2 * ll + attr(object$penalty, "k") * object$df

    ## order
    pi <- order(aic)
    pi <- array(c((pi - 1) %%  object$nbest + 1,
                  (pi - 1) %/% object$nbest + 1),
                dim = c(length(pi), 2))
    nok <- is.na(aic[pi])
    pi <- pi[!nok          , , drop = FALSE]
    pi <- pi[1:object$nbest, , drop = FALSE]

    x.names <- variable.names(object, .full = TRUE)
    .cnames <- paste(1:object$nbest, ".", sep = "")
    ## df
    object$df <- array(object$df[pi], dim = object$nbest,
                       dimnames = list(.cnames))
    ## rss
    object$rss <- array(object$rss[pi], dim = object$nbest,
                        dimnames = list(.cnames))
    ## val
    object$val <- array(aic[pi], dim = object$nbest,
                        dimnames = list(.cnames))
    ## which table
    sel <- cbind(rep(1:object$nvar, times = nrow(pi)   ),
                 rep(pi[, 1]      , each  = object$nvar),
                 rep(pi[, 2]      , each  = object$nvar))
    object$which <- array(object$which[sel], dim = c(object$nvar, object$nbest),
                          dimnames = list(x.names, .cnames))

    ## class
    class(object) <- "lmSelect"
    
    ## done
    object
}


## workhorse function
##
## Args:
##   x         - (matrix) model matrix
##   y         - (numeric[]) response variable
##   weights   - (numeric[]) 
##   offset    - (numeric[])
##   include   - (integer[]) regressors to force in
##   exclude   - (integer[]) regressors to force out
##   penalty   - ("AIC"|"BIC"|numeric) penalty per parameter
##   tolerance - (numeric) tolerance
##   pradius   - (integer) preordering radius
##   nbest     - (integer) number of best subsets
##   ...       - ignored
##   .algo     - (character) for interal use
##
## Rval:  (lmSelect)
##   call      - (call)
##   nobs      - (integer)
##   nvar      - (integer)
##   weights   - (numeric[])
##   offset    - (numeric[])
##   intercept - (logical)
##   include   - (integer[])
##   exclude   - (integer[])
##   nmin      - (integer)
##   nmax      - (integer)
##   penalty   - (numeric)
##   tolerance - (numeric[])
##   nbest     - (integer)
##   df        - (integer[])
##   rss       - (array)
##   val       - (array)
##   which     - (array)
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
##   'xbba*':  experimental PBBA
##
lmSelect.fit <- function (x, y, weights = NULL, offset = NULL,
                          include = NULL, exclude = NULL, penalty = "BIC",
                          tolerance = 0, pradius = NULL, nbest = 1, ...,
                          .algo = "hpbba")
{
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
    nmin <- mark + 1
    nmax <- size

    ## preordering radius
    ## TODO: validate
    if (is.null(pradius)) {
        pmin <- 8
    } else {
        pmin <- size - mark - pradius + 1
    }

    ## penalty
    if (tolower(penalty) == "aic") {
        penalty <- structure("AIC", k = 2)
    } else if (tolower(penalty) == "bic") {
        penalty <- structure("BIC", k = log(nobs))
    } else if (is.numeric(penalty)) {
        penalty <- structure("AIC", k = penalty)
    } else {
        stop ("invalid 'penalty'")
    }


    ## tolerance
    ## TODO: validate

    ## call C
    C_v <- which - 1
    C_xy <- cbind(x[, which], y)
    C_tau <- tolerance + 1

    C_index <- rep(0, nbest)
    C_df <- rep(-1, nbest)
    C_rss <- rep(0, nbest)
    C_val <- rep(0, nbest)
    C_which <- rep(0, nbest * nvar)

    C_penalty <- attr(penalty, "k")

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
        val     = as.double(C_val),
        which   = as.logical(C_which),
        penalty = as.double(C_penalty),
        tau     = as.double(C_tau),
        info    = integer(1),
        nodes   = integer(1)
    )

    C_rval <- do.call(".C", c("R_lmSelect", C_args))
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
                 penalty   = penalty,
                 tolerance = NULL,
                 nbest     = nbest,
                 df        = NULL,
                 rss       = NULL,
                 val       = NULL,
                 which     = NULL,
                 .info     = C_rval$info,
                 .nodes    = C_rval$nodes)
    class(rval) <- "lmSelect"

    ## extract value and subsets
    .cnames <- format.ordinal(1:nbest)
    ## tolerance
    rval$tolerance <- tolerance
    ## df
    rval$df <- array(C_rval$df, dim = nbest, dimnames = list(.cnames))
    not.valid <- rval$df < 0
    rval$df[not.valid] <- NA
    ## rss
    rval$rss <- array(C_rval$rss, dim = nbest, dimnames = list(.cnames))
    rval$rss[not.valid] <- NA
    ## val
    S <- sum(log(w))
    rval$val <- array(C_rval$val, dim = nbest, dimnames = list(.cnames)) - S
    rval$val[not.valid] <- NA
    ## which
    rval$which <- array(C_rval$which, dim = c(nvar, nbest),
                        dimnames = list(x.names, .cnames))
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
##   object    - (formula) an object that can be coerced to class 'formula'
##   data      - (data.frame)
##   subset    - (vector) subset of observations
##   weights   - (numeric[])
##   na.action - (function)
##   model     - (logical)
##   x         - (logical)
##   y         - (logical)
##   contrasts - (list)
##   offset    - (numeric[])
##   ...       - forwarded to 'lmSelect.fit'
##
## Rval:  (lmSelect) see 'lmSelect.fit'
##
lmSelect <- function (formula, data, subset, weights, na.action,
                      model = TRUE, x = FALSE, y = FALSE,
                      contrasts = NULL, offset, ...) {
    ## coerce 'lmSubsets' object
    if (inherits(formula, "lmSubsets")) {
        return (lmSubsets.select(formula, ...))
    }

    ## construct 'lmSelect' object
    ret.m <- model;  model <- NULL
    ret.x <- x;  x <- NULL
    ret.y <- y;  y <- NULL

    ## keep call (of generic)
    cl <- match.call()
    cl[[1]] <- as.name("lmSelect")

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
        stop("'weights' must be a numeric vector")
    }

    ## offset
    offset <- as.vector(model.offset(mf))
    if(!is.null(offset) && (length(offset) != NROW(y))) {
        stop("length of 'offset' must equal number of observations")
    }

    ## fit subsets
    rval <- lmSelect.fit(x, y, w, offset = offset, ...)

    ## return value
    class(rval) <- "lmSelect"
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


## print method for 'lmSelect' objects
##
## Args:
##   x      - (lmSelect)
##   ...    - ignored
##
## Rval:  (lmSelect) invisible
##
print.lmSelect <- function (x, ...)
{
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
                       format(x$nmin),
                       format(x$nmax),
                       paste(x$penalty, " (k = ", format(attr(x$penalty, "k"), nsmall = 2), ")"),
                       format(x$tolerance),
                       format(x$nbest)),
                     ncol = 1)
    colnames(val) <- ""
    rownames(val) <- paste("  ",
                           c("N observations", "N regressors", "Weights",
                             "Offset", "Intercept", "N min", "N max", "Value",
                             "Tolerance", "N best"),
                           ":")
    print(val, quote = FALSE)


    ## fit
    catln()
    catln("Model fit:")
    fit <- format(rbind(x$df, x$rss, x$val), nsmall = 2)
    rownames(fit) <- paste("  ", c("df", "Deviance", "Value"))
    print(fit, quote = FALSE)


    ## which
    catln()
    catln("Which:")
    which <- apply(x$which[, , drop = FALSE], 1, format.which)
    which <- matrix(which, dimnames = list(names(which), "best"))
    rownames(which)[x$include] <- paste("+", rownames(which)[x$include])
    rownames(which)[x$exclude] <- paste("-", rownames(which)[x$exclude])
    rownames(which) <- paste("  ", rownames(which))
    print(which, quote = FALSE)


    catln()

    ## done
    invisible(x)
}


## plot method for 'lmSelect' objects
##
## Args:
##   x      - (lmSelect)
##   ...    - forwarded to plotting functions
##   xlim   - (numeric[])
##   ylim1  - (numeric[])
##   ylim2  - (numeric[])
##   legend - (character[])
##
## Rval:  (lmSelect) invisible
##
## Graphical arguments are passed to 'plot.default'.
##
plot.lmSelect <- function (x, ..., xlim = NULL, ylim1 = NULL, ylim2 = NULL,
                           legend) {
    localPlot <- function (object, main, sub = NULL, xlab, ylab, type = "o",
                           lty = c(3, 1), pch = c(21, 16), col = "black",
                           bg = "white", ...) {
        if (missing(main)) main <- "Best subsets"
        if (missing(xlab)) xlab <- "Best subset"
        if (missing(ylab)) ylab <- c("Deviance", "Value")

        type <- rep(type, length = 2)
        lty  <- rep(lty , length = 2)
        pch  <- rep(pch , length = 2)
        col  <- rep(col , length = 2)
        bg   <- rep(bg  , length = 2)

        x <- seq(object$nbest)
        y1 <- object$rss
        y2 <- object$val

        if (is.null(xlim)) xlim <- range(x)
        if (is.null(ylim1)) ylim1 <- range(y1[is.finite(y1)])
        if (is.null(ylim2)) ylim2 <- range(y2[is.finite(y2)])

        par(mar = c(5, 4, 4, 4) + 0.1)

        plot(x, y1, main = main, sub = sub, xlab = xlab, ylab = ylab[1],
             type = type[1], xlim = xlim, ylim = ylim1, lty = lty[1],
             pch = pch[1], col = col[1], bg = bg[1], ...)

        if (!is.null(legend)) {
            legend("topleft", legend = legend, lty = lty, pch = pch, col = col,
                   pt.bg = bg, bty = "n")
        }

        plot.window(xlim = xlim, ylim = ylim2)

        axis(side = 4, at = pretty(ylim2))
        mtext(ylab[2], side = 4, line = 3)

        lines(x, y2, type = type[1], lty = lty[2], pch = pch[2], col = col[2],
              bg = bg[2], ...)
    }

    if (inherits(x, "summary.lmSelect")) {
        if (missing(legend)) legend <- c("Deviance (RSS)", paste("Value (", x$penalty, ")", sep = ""))

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
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (character[])
##   variable names
##
variable.names.lmSelect <- function (object, best = 1, ...) {
    ## full model
    x.names <- dimnames(object$which)[[1]]

    ## subset model
    x.names[object$which[, best]]
}


## extract formula
##
## Args:
##   x    - (lmSelect)
##   best - (integer)
##   ...  - ignored
##
## Rval:  (formula)
##
formula.lmSelect <- function (x, best, ...) {
    if (missing(best)) {
        f <- formula(x$terms)

        return (f)
    }

    e <- new.env();
    e$x <- model.matrix(x, best = best)
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
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (character[])
##   variable names
##
model.frame.lmSelect <- function(formula, best, ...)
{
    mf <- formula[["model"]]

    if (!missing(best) || is.null(mf)) {
        cl <- formula$call
        m <- c("formula", "data", "subset",
               "weights", "na.action", "offset")
        m <- match(m, names(cl), 0L)
        cl <- cl[c(1L, m)]
        cl$drop.unused.levels <- TRUE
        cl[[1L]] <- quote(stats::model.frame)
        cl$xlev <- formula$xlevels
        cl$formula <- stats::formula(formula, best = best)

        env <- environment(formula$terms)
        if (is.null(env)) {
            env <- parent.frame()
        }

        mf <- eval(cl, env)
    }

    mf
}


## extract model response
##
## Args:
##   x    - (lmSubsets)
##   ...  - ignored
##
## Rval:  (formula)
##
model.response.lmSelect <- function (data, ...) {
    y <- data[["y"]]
    if (is.null(y)) {
        mf <- model.frame(data)
        y <- model.response(mf)
    }

    y
}


## extract model matrix
##
## Args:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (matrix)
##
model.matrix.lmSelect <- function (object, best, ...) {
    ## full model
    x <- object[["x"]]
    if (is.null(x)) {
        data <- model.frame(object)
        x <- model.matrix.default(object, data = data, contrasts.arg = object$contrasts)
    }

    if (missing(best)) {
        return (x)
    }

    ## subset model
    x[, object$which[, best]]
}


## refit
##
## Args:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - forwarded to 'lm'
##
## Rval:  (lm)
##
refit.lmSelect <- function (object, best = 1, ...) {
    f <- formula(object, best = best)

    lm(f, ...)
}


## extract ceofficients
##
## Args:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   This method refits the 'lmSelect' object and calls 'coef' on the
##   obtained 'lm' object.  In some circumstances it might be more
##   efficient to call 'refit' directly and execute 'coef' on the
##   obtained 'lm' object.
##
coef.lmSelect <- function (object, best = 1, ...) {
    coef(refit(object, best = best, model = FALSE))
}


## extract variance-covariance matrix
##
## Args:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (matrix)
##
## Note: 'refit'
##   See 'coef.lmSelect'.
##
vcov.lmSelect <- function (object, best = 1, ...) {
    vcov(refit(object, best = best, model = FALSE))
}


## extract fitted values
##
## Args:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.lmSelect'.
##
fitted.lmSelect <- function (object, best = 1, ...) {
  fitted(refit(object, best = best, model = FALSE))
}


## extract residuals
##
## Args:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.lmSelect'.
##
residuals.lmSelect <- function (object, best = 1, ...) {
  residuals(refit(object, best = best, model = FALSE))
}


## extract deviance
##
## Args:
##   object        - (lmSelect)
##   best          - (integer[])
##   ...           - ignored
##
## Rval:  (numeric[])
##
## Returns the RSS for the specified subset(s).
##
deviance.lmSelect <- function (object, best = 1, ...) {
    ## extract RSS
    object$rss[best]
}


## extract log-likelihood (base method)
##
## Args:
##   object - (lmSelect)
##   best   - (integer[])
##   ...    - ignored
##
## Rval: (logLik)
##
## Note:
##   Can handle multiple objects.
##
logLik.lmSelect <- function (object, best = 1, ...) {
    ## weights
    N <- object$nobs
    if (is.null(w <- object$weights)) {
        S <- 0
    } else {
        S <- sum(log(w[w != 0]))
    }

    ## extract rss
    rss <- deviance(object, best = best)

    ## done
    structure(0.5 * (S - N * (log(2 * pi) + 1 - log(N) + log(rss))),
              df       = object$df[best]                           ,
              nobs     = N                                         ,
              class    = "logLik"                                  )
}


## compute AIC
##
## Args:
##   object - (lmSelect)
##   best   - (integer[])
##   ...    - ignored
##   k      - (numeric)
##
## Rval: (numeric|data.frame)
##
AIC.lmSelect <- function (object, best = 1, ..., k = 2) {
    ## extract log-likelihoods
    ll <- logLik(object, best = best)
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
##   object - (lmSelect)
##   best   - (integer[])
##   ...    - ignored
##
## Rval: (numeric|data.frame)
##
BIC.lmSelect <- function (object, best = 1, ...) {
    ## extract log-likelihoods
    ll <- logLik(object, best = best)
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
