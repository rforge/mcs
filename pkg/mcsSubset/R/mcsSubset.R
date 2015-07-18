##
## File:  mcsSubset.R
##



################
## GENERATORS ##
################


## standard formula interface
##
## Args:
##   formula   - (formula)
##   ...       - forwarded to 'lm' and 'mcsSubset.lm'
##   lm        - (logical) if 'true', compute 'lm' component
##
## Rval:  (mcsSubset)
##   See 'mcsSubset.default'.
##
## NOTE:  'lm'
##   If 'lm' is 'FALSE', the returned 'lm' component is
##   an "empty" mockup.
##
mcsSubset.formula <- function (formula, ..., lm = FALSE) {
    ret.lm <- lm;  lm <- NULL

    ## keep call (of generic)
    call <- match.call()
    call[[1]] <- as.name("mcsSubset")

    ## lm call
    lm.call <- match.call(expand.dots = TRUE)
    m <- match(c("include", "exclude", "size", "penalty", "tolerance",
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
        m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
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
    cl[[1L]] <- as.name("mcsSubset")
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
##   object       - (lm)
##   ...          - forwarded to 'mcsSubset.default'
##
## Rval:  (mcsSubset)
##   See 'mcsSubset.default'.
##
mcsSubset.lm <- function (object, ...) {
    ## keep call (of generic)
    call <- match.call()
    call[[1]] <- as.name("mcsSubset")

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
    if (is.null(w <- model.weights(mf))) w <- rep(1, NROW(x))
    if (is.null(o <- model.offset(mf))) o <- rep(0, NROW(x))
    wok <- w != 0;  wnz <- w[wok]
    x <- sqrt(wnz) * x[wok, , drop = FALSE]
    y <- sqrt(wnz) * y[wok,   drop = FALSE] - o[wok]

    ## forward call
    rval <- mcsSubset(x, y, ...)

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
##   include   - (integer[]) regressors to force in
##   exclude   - (integer[]) regressors to force out
##   size      - (integer[]) subset sizes
##   penalty   - (numeric) AIC penalty
##   tolerance - (numeric[]) tolerance
##   pradius   - (integer) preordering radius
##   nbest     - (integer) number of best subsets
##   ...       - ignored
##   .algo     - (character) for interal use
##
## Rval:  (mcsSubset)
##   call      - (call)
##   nobs      - (integer)
##   nvar      - (integer)
##   include   - (integer[])
##   exclude   - (integer[])
##   size      - (integer[])
##   intercept - (logical)
##   penalty   - (numeric)
##   tolerance - (numeric[])
##   nbest     - (integer)
##   rss       - (array)
##   aic       - (array)
##   which     - (array)
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
mcsSubset.default <- function (object, y, include = NULL, exclude = NULL,
                               size = NULL, penalty = 0, tolerance = 0,
                               pradius = NULL, nbest = 1, ..., .algo = "hbba")
{
    ## keep call (of generic)
    call <- match.call()
    call[[1]] <- as.name("mcsSubset")

    ## model matrix
    x <- as.matrix(object);  object <- NULL

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
    intercept <- isTRUE(all.equal(as.vector(x[, 1]), rep(1, nobs)))
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
    if (penalty == 0) {
        tolerance <- rep(tolerance, length.out = nvar)
        tolerance[-size] <- +Inf
    } else {
        tolerance <- tolerance[1]
    }

    ## call C
    C_mark <- size.min - 1
    C_size <- size.max
    C_v <- which - 1
    C_xy <- cbind(x[, which], y)
    C_tau <- tolerance + 1
    if (penalty == 0) {
        C_index <- rep(0, nvar * nbest)
        C_rss <- rep(0, nvar * nbest)
        C_aic <- double(0)
        C_which <- rep(0, nvar * nbest * nvar)

        C_tau <- C_tau[size]
        C_tau[C_tau == +Inf] <- .Machine$double.xmax
        C_tau <- c(0, C_tau)
    } else {
        C_index <- rep(0, nbest)
        C_rss <- rep(0, nbest)
        C_aic <- rep(0, nbest)
        C_which <- rep(0, nbest * nvar)
    }

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
        aic     = as.double(C_aic),
        which   = as.logical(C_which),
        penalty = as.double(penalty),
        tau     = as.double(C_tau),
        nodes   = integer(1),
        info    = integer(1)
    )

    C_rval <- do.call(".C", c("R_mcsSubset", C_args))

    ## return value
    rval <- list(call      = call,
                 nobs      = nobs,
                 nvar      = nvar,
                 include   = include,
                 exclude   = exclude,
                 size      = NULL,
                 intercept = intercept,
                 penalty   = penalty,
                 tolerance = NULL,
                 nbest     = nbest,
                 rss       = NULL,
                 aic       = NULL,
                 which     = NULL,
                 .nodes = C_rval$nodes)
    class(rval) <- "mcsSubset"

    ## extract value & subsets
    if (penalty == 0) {
        .cnames <- paste(1:nbest, ".", sep = "")
        ## size
        rval$size <- size
        ## tolerance
        rval$tolerance <- array(tolerance, dim = nvar,
                                dimnames = list(1:nvar))
        ## rss
        rval$rss <- array(C_rval$rss, dim = c(nbest, nvar),
                          dimnames = list(.cnames, 1:nvar))
        rval$rss[, -rval$size] <- NA
        rval$rss[rval$rss == .Machine$double.xmax] <- NA
        ## which
        rval$which <- array(C_rval$which, dim = c(nvar, nbest, nvar),
                            dimnames = list(x.names, .cnames, 1:nvar))
        rval$which[, , -rval$size] <- NA
    } else {
        .cnames <- paste(1:nbest, ".", sep = "")
        ## tolerance
        rval$tolerance <- tolerance
        ## rss
        rval$rss <- array(C_rval$rss, dim = nbest, dimnames = list(.cnames))
        ## aic
        rval$aic <- array(C_rval$aic, dim = nbest, dimnames = list(.cnames))
        ## which
        rval$which <- array(C_rval$which, dim = c(nvar, nbest),
                            dimnames = list(x.names, .cnames))
        ## size
        rval$size <- array(apply(rval$which, 2, sum), dim = nbest,
                           dimnames = list(.cnames))
    }

    ## done
    rval
}


## print method for 'mcsSubset' objects
##
## Args:
##   x      - (mcsSubset)
##   digits - (integer)
##   ...    - ignored
##
## Rval:  (mcsSubset) invisible
##
print.mcsSubset <- function (x, digits = NULL, ...)
{
    catln <- function (...) base::cat(..., "\n", sep = "")

    if (is.null(digits)) {
        digits <- max(3, getOption("digits") - 3)
    }

    catln()
    catln("Call:")
    catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))

    val <- as.matrix(c(format(x$nvar),
                       if (x$intercept) "Yes" else "No",
                       paste(x$include, collapse = " "),
                       paste(x$exclude, collapse = " "),
                       paste(x$size, collapse = " "),
                       if (x$penalty == 0) {
                           "RSS"
                       } else {
                           paste("AIC (k = ", format(x$penalty, digits = digits),
                                 ")", sep = "")
                       },
                       paste(format(x$tolerance, digits = digits, trim = TRUE),
                             collapse = " "),
                       format(x$nbest)),
                     ncol = 1)
    colnames(val) <- ""
    rownames(val) <- c("Total regressors:", "Intercept:", "Include:", "Exclude:",
                       "Size:", "Criterion:", "Tolerance:", "N best:")

    print(val, quote = FALSE)

    invisible(x)
}


## plot method for 'mcsSubset' objects
##
## Args:
##   x      - (mcsSubset)
##   type   - (character)
##   main   - (character)
##   xlab   - (character)
##   ylab   - (character)
##   col    - (integer[]|character[])
##   lty    - (integer)
##   ...    - passed to 'deviance.mcsSubset' or 'AIC.mcsSubset'
##   legend - (logical)
##
## Rval:  (mcsSubset) invisible
##
## Graphical arguments are passed to 'plot.default'.
##
plot.mcsSubset <- function (x, type = "b", main = "Deviance",
                            xlab = NULL, ylab = "", col = "blue",
                            lty = 1, ..., legend = TRUE) {
    if (x$penalty == 0) {
        ## rss
        rss <- deviance(x, ...)
        ## sub title
        sub <- "RSS"
        ## xlab
        if (is.null(xlab)) {
            xlab <- "Number of regressors in model"
        }
        ## plot
        plot(names(rss), rss, type = type, main = main, sub = sub,
             xlab = xlab, ylab = ylab, col = col, lty = lty)
        ## legend
        if (legend) {
            legend("topright", "RSS", col = col, lty = lty, bty = "n")
        }
    } else {
        digits <- max(3, getOption("digits") - 3)
        ## aic
        aic <- AIC(x, ...)$AIC
        ## sub title
        sub <- paste("AIC (k = ", format(x$penalty, digits = digits), ")", sep = "")
        ## xlab
        if (is.null(xlab)) {
            xlab <- "Best"
        }
        ## best
        best <- 1:x$nbest
        ## plot
        plot(best, aic, type = type, main = main, sub = sub,
             xlab = xlab, ylab = ylab, col = col, lty = lty)
        ## labels (subset size)
        text(best, aic, labels = paste0("(", best, ")"), ## FIXME: was "size"
             adj = c(0, 1.5), cex = 0.8)
        ## legend
        if (legend) {
            legend("topleft", "AIC (size)", lty = lty, col = col, bty = "n")
        }
    }

    ## done
    invisible(x)
}



#############
## METHODS ##
#############


## extract variable names
##
## Args:
##   object - (mcsSubset)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##   .full  - (logical) for internal use
##   .neg   - (logical) for internal use
##
## Rval:  (character[])
##   variable names
##
variable.names.mcsSubset <- function (object, size = NULL, best = 1, ...,
                                      .full = FALSE, .neg = FALSE) {
    if (is.null(object$.lm)) {
        stop ("'mcsSubset' object does not have a 'lm' component")
    }

    ## full model
    x.names <- variable.names(object$.lm)
    if (.full) {
        return (x.names)
    }
    ## which
    if (object$penalty == 0) {
        wi <- object$which[, best, size]
    } else {
        wi <- object$which[, best]
    }
    ## negate
    wi <- xor(wi, .neg)
    ## select
    x.names <- x.names[wi]

    ## done
    x.names
}


## extract formula
##
## Args:
##   x   - (mcsSubset)
##   ... - forwarded to 'variable.names.mcsSubset'
##
## Rval:  (formula)
##
formula.mcsSubset <- function (x, ...) {
    if (is.null(x$.lm)) {
        stop ("'mcsSubset' object does not have a 'lm' component")
    }

    ## full model
    f <- formula(x$.lm)
    ## variable names
    x.names <- variable.names(x, ..., .neg = TRUE)
    ## update formula
    f <- update(f, paste(".~.", paste(c("", x.names), collapse = "-")))

    ## done
    f
}


## extract model frame
##
## Args:
##   formula       - (mcsSubset)
##   ...           - forwarded to 'formula.mcsSubset'
##
## Rval:  (data.frame)
##
model.frame.mcsSubset <- function (formula, ...) {
    if (is.null(formula$.lm)) {
        stop ("'mcsSubset' object does not have a 'lm' component")
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
##   object        - (mcsSubset)
##   size          - (integer)
##   best          - (integer)
##   ...           - ignored
##
## Rval:  (matrix)
##
model.matrix.mcsSubset <- function (object, size = NULL, best = 1, ...) {
    if (is.null(object$.lm)) {
        stop ("'mcsSubset' object does not have a 'lm' component")
    }

    ## full model
    x <- model.matrix(object$.lm)
    ## which
    if (object$penalty == 0) {
        wi <- object$which[, best, size]
    } else {
        wi <- object$which[, best]
    }
    ## select
    x <- x[, wi]

    ## done
    x
}


## refit
##
## Args:
##   object - (mcsSubset)
##   ...    - forwarded to 'formula.mcsSubset'
##
## Rval:  (lm)
##
refit.mcsSubset <- function (object, ...) {
    if (is.null(object$.lm)) {
        stop ("'mcsSubset' object does not have a 'lm' component")
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
##   object - (mcsSubset)
##   ...    - forwarded to 'refit.mcsSubset'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   This method refits the 'mcsSubset' object and calls 'coef' on the
##   obtained 'lm' object.  In some circumstances it might be more
##   efficient to call 'refit' directly and execute 'coef' on the
##   obtained 'lm' object.
##
coef.mcsSubset <- function (object, ...) {
    coef(refit(object, ...))
}


## extract variance-covariance matrix
##
## Args:
##   object - (mcsSubset)
##   ...    - forwarded to 'refit.mcsSubset'
##
## Rval:  (matrix)
##
## Note: 'refit'
##   See 'coef.mcsSubset'.
##
vcov.mcsSubset <- function (object, ...) {
    vcov(refit(object, ...))
}


## extract fitted values
##
## Args:
##   object - (mcsSubset)
##   ...    - forwarded to 'refit.mcsSubset'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.mcsSubset'.
##
fitted.mcsSubset <- function (object, ...) {
  fitted(refit(object, ...))
}


## extract residuals
##
## Args:
##   object - (mcsSubset)
##   ...    - forwarded to 'refit.mcsSubset'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.mcsSubset'.
##
residuals.mcsSubset <- function (object, ...) {
  residuals(refit(object, ...))
}


## extract deviance
##
## Args:
##   object        - (mcsSubset)
##   size          - (integer[])
##   best          - (integer[])
##   ...           - ignored
##
## Rval:  (numeric[])
##
## Returns the RSS for the specified subset(s).
##
deviance.mcsSubset <- function (object, size = NULL, best = 1, ...) {
    ## compute indices
    if (object$penalty == 0) {
        if (is.null(size)) size <- object$size
        d <- object$rss[best, size]
        names(d) <- size
    } else {
        d <- object$rss[best]
    }

    ## done
    d
}


## extract log-likelihood (base method)
##
## Args:
##   object - (mcsSubset)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##   df     - (integer[]) degrees of freedom
##
## Rval: (logLik)
##
## Note:
##   Can handle multiple objects.
##
logLik.mcsSubset <- function (object, size = NULL, best = 1, ..., df) {
    if (is.null(object$.lm)) {
        stop ("'mcsSubset' object does not have a 'lm' component")
    }

    ## degrees of freedom
    if (object$penalty == 0) {
        if (is.null(size)) size <- object$size
        df <- rep(size + 1, each = length(best))
    } else {
        df <- object$size[best] + 1
    }

    ## weights
    if(is.null(w <- weights(object$.lm))) {
        nobs <- object$nobs
        sw <- 0
    } else {
        wnz <- w[w > 0]
        nobs <- sum(wnz)
        sw <- sum(log(wnz))
    }

    ## extract rss
    rss <- deviance(object, size = size, best = best)

    ## done
    structure(0.5 * (sw - nobs * (log(2 * pi) + 1 - log(nobs) + log(rss))),
              df = df, nobs = nobs, class = "logLik")
}


## compute AIC
##
## Args:
##   object - (mcsSubset)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##   k      - (integer)
##
## Rval: (numeric|data.frame)
##
AIC.mcsSubset <- function (object, size = NULL, best = 1, ..., k = NULL) {
    if (object$penalty == 0) {
        if (is.null(size)) size <- object$size
        if (is.null(k)) k <- 2
    } else {
        if (is.null(k)) k <- object$penalty
    }

    ## extract log-likelihoods
    ll <- logLik(object, size = size, best = best)
    ## compute AICs
    aic <- AIC(ll, k = k)
    ## data frame?
    if (length(aic) > 1) {
        aic <- data.frame(df = attr(ll, "df"), AIC = aic)
    }

    ## done
    aic
}
