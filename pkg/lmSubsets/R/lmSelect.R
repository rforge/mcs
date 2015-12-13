##
## File:  lmSelect.R
##



################
## GENERATORS ##
################


## construct 'lmSelect' from 'lmSubsets' object
##
## Args:
##   object   - (lmSubsets)
##   ...      - ignored
##   penalty  - ("AIC"|"BIC"|numeric) penalty per parameter
##
## Rval:  (lmSelect)
##   See 'lmSelect.default'.
##
lmSelect.lmSubsets <- function (object, ..., penalty = "AIC") {
    paste <- function (..., sep = "") base::paste(..., sep = sep)

    ## tolerance
    tolerance <- object$tolerance[object$size]
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
        S <- sum(log(w))
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
    .cnames <- paste(1:object$nbest, ".")
    ## df
    object$df <- array(object$df[pi], dim = object$nbest,
                       dimnames = list(.cnames))
    ## rss
    object$rss <- array(object$rss[pi], dim = object$nbest,
                        dimnames = list(.cnames))
    ## aic
    object$aic <- array(aic[pi], dim = object$nbest,
                        dimnames = list(.cnames))
    ## which table
    sel <- cbind(rep(1:object$nvar, times = nrow(pi)   ),
                 rep(pi[, 1]      , each  = object$nvar),
                 rep(pi[, 2]      , each  = object$nvar))
    object$which <- array(object$which[sel], dim = c(object$nvar, object$nbest),
                          dimnames = list(x.names, .cnames))
    ## size
    object$size <- NULL

    ## class
    class(object) <- "lmSelect"
    
    ## done
    object
}


## standard formula interface
##
## Args:
##   formula   - (formula)
##   ...       - forwarded to 'lm' and 'lmSelect.lm'
##   lm        - (logical) if 'true', compute 'lm' component
##
## Rval:  (lmSelect)
##   See 'lmSelect.default'.
##
## NOTE:  'lm'
##   If 'lm' is 'FALSE', the returned 'lm' component is
##   an "empty" mockup.
##
lmSelect.formula <- function (formula, ..., lm = FALSE) {
    ret.lm <- lm;  lm <- NULL

    ## keep call (of generic)
    call <- match.call()
    call[[1]] <- as.name("lmSelect")

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
    cl[[1L]] <- as.name("lmSelect")
    cl$formula <- cl$y <- cl$lm <- NULL
    cl$object <- object
    rval <- eval(cl, parent.frame())

    ## return value
    rval$call <- call
    if (!ret.lm) rval$lm <- NULL
    if (!ret.m) rval$.lm$model <- NULL
    if (!ret.x) rval$.lm$x <- NULL
    if (!ret.y) rval$.lm$y <- NULL

    ## done
    rval
}


## interface for fitted lm regressions
##
## Args:
##   object       - (lm)
##   ...          - forwarded to 'lmSelect.default'
##   penalty      - ("AIC"|"BIC"|numeric) penalty per parameter
##
## Rval:  (lmSelect)
##   See 'lmSelect.default'.
##
lmSelect.lm <- function (object, ..., penalty = "AIC") {
    ## keep call (of generic)
    call <- match.call()
    call[[1]] <- as.name("lmSelect")

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
    rval <- lmSelect(x, y, weights = w, offset = o, penalty = penalty, ...)

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
##   include   - (integer[])
##   exclude   - (integer[])
##   size      - (integer[])
##   intercept - (logical)
##   penalty   - (numeric)
##   tolerance - (numeric[])
##   nbest     - (integer)
##   df        - (integer[])
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
lmSelect.default <- function (object, y, weights = NULL, offset = NULL,
                              include = NULL, exclude = NULL, penalty = "AIC",
                              tolerance = 0, pradius = NULL, nbest = 1, ...,
                              .algo = "hbba")
{
    ## keep call (of generic)
    call <- match.call()
    call[[1]] <- as.name("lmSelect")

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
    size <- size.min:size.max

    ## preordering radius
    ## TODO: validate
    if (is.null(pradius)) {
        pmin <- 8
    } else {
        pmin <- size.max - pradius
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
    C_mark <- size.min - 1
    C_size <- size.max
    C_v <- which - 1
    C_xy <- cbind(x[, which], y)
    C_tau <- tolerance + 1

    C_index <- rep(0, nbest)
    C_rss <- rep(0, nbest)
    C_aic <- rep(0, nbest)
    C_which <- rep(0, nbest * nvar)

    C_penalty <- attr(penalty, "k")

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
        penalty = as.double(C_penalty),
        tau     = as.double(C_tau),
        nodes   = integer(1),
        info    = integer(1)
    )

    C_rval <- do.call(".C", c("R_lmSelect", C_args))

    ## return value
    rval <- list(call      = call,
                 nobs      = nobs,
                 nvar      = nvar,
                 weights   = weights[wok],
                 offset    = offset[wok],
                 include   = include,
                 exclude   = exclude,
                 size      = size,
                 intercept = intercept,
                 penalty   = penalty,
                 tolerance = NULL,
                 nbest     = nbest,
                 df        = NULL,
                 rss       = NULL,
                 aic       = NULL,
                 which     = NULL,
                 .nodes = C_rval$nodes)
    class(rval) <- "lmSelect"

    ## extract value & subsets
    .cnames <- paste(1:nbest, ".", sep = "")
    ## tolerance
    rval$tolerance <- tolerance
    ## rss
    rval$rss <- array(C_rval$rss, dim = nbest, dimnames = list(.cnames))
    ## aic
    S <- sum(log(w))
    rval$aic <- array(C_rval$aic, dim = nbest, dimnames = list(.cnames)) - S
    ## which
    rval$which <- array(C_rval$which, dim = c(nvar, nbest),
                        dimnames = list(x.names, .cnames))
    ## df
    rval$df <- array(colSums(rval$which) + 1, dim = nbest,
                     dimnames = list(.cnames))

    ## done
    rval
}


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
                       paste(x$size, collapse = " "),
                       paste(x$penalty, " (k = ", format(attr(x$penalty, "k"), nsmall = 2), ")"),
                       format(x$tolerance, nsmall = 2),
                       format(x$nbest)),
                     ncol = 1)
    colnames(val) <- ""
    rownames(val) <- paste("  ",
                           c("Number of observations", "Number of regressors",
                             "Weights", "Offset", "Intercept", "Include", "Exclude",
                             "Subset sizes", "Value", "Tolerance", "N best"),
                           ":")
    print(val, quote = FALSE)


    ## fit
    catln()
    catln("Model fit:")
    fit <- format(rbind(x$df, x$rss, x$aic), nsmall = 2)
    rownames(fit) <- paste("  ", c("df", "Deviance", "Value"))
    print(fit, quote = FALSE)
    catln()


    ## done
    invisible(x)
}


## plot method for 'lmSelect' objects
##
## Args:
##   x      - (lmSelect)
##   ...    - ignored
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
    localPlot <- function (object, main, sub = NULL, xlab, ylab,
                           type, lty, pch, col, bg, ...) {
        if (missing(main)) main <- "Best subsets"
        if (missing(xlab)) xlab <- "Best subset"
        if (missing(ylab)) ylab <- c("Deviance", "Value")

        type <- if (missing(type)) c("o", "o")         else rep(type, length = 2)
        lty  <- if (missing(lty) ) c(3, 1)             else rep(lty , length = 2)
        pch  <- if (missing(pch) ) c(21, 21)           else rep(pch , length = 2)
        col  <- if (missing(col) ) c("black", "red")   else rep(col , length = 2)
        bg   <- if (missing(bg)  ) c("white", "white") else rep(bg  , length = 2)
        
        x <- seq(object$nbest)
        y1 <- object$rss
        y2 <- object$aic

        if (is.null(xlim)) xlim <- range(x)
        if (is.null(ylim1)) ylim1 <- range(y1[is.finite(y1)])
        if (is.null(ylim2)) ylim2 <- range(y2[is.finite(y2)])

        par(mar = c(5, 4, 4, 4) + 0.1)

        plot(x, y1, main = main, sub = sub, xlab = xlab, ylab = ylab[1],
             type = type[1], xlim = xlim, ylim = ylim1, lty = lty[1],
             pch = pch[1], col = col[1], bg = bg[1], ...)

        if (!is.null(legend)) {
            legend("topleft", legend = legend, lty = lty,
                   pch = pch, col = col, pt.bg = bg, bty = "n")
        }

        plot.window(xlim = xlim, ylim = ylim2)

        axis(side = 4, at = pretty(ylim2))
        mtext(ylab[2], side = 4, line = 3)

        lines(x, y2, type = type[2], lty = lty[2], pch = pch[2],
              col = col[2], bg = bg[2], ...)
    }

    if (missing(legend)) legend <- c("Deviance (RSS)", paste("Value (", x$penalty, ")", sep = ""))

    localPlot(x, ...)

    invisible(x)
}



#############
## METHODS ##
#############


## extract variable names
##
## Args:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##   .full  - (logical) for internal use
##   .cmpl  - (logical) for internal use
##
## Rval:  (character[])
##   variable names
##
variable.names.lmSelect <- function (object, best = 1, ...,
                                     .full = FALSE, .cmpl = FALSE) {
    ## lm
    if (is.null(object$.lm)) {
        stop ("'lmSelect' object does not have a 'lm' component")
    }

    ## full model
    x.names <- variable.names(object$.lm)
    if (.full) {
        return (x.names)
    }

    ## which
    wi <- object$which[, best]
    ## complementary model
    wi <- xor(wi, .cmpl)
    ## select
    x.names <- x.names[wi]

    ## done
    x.names
}


## extract formula
##
## Args:
##   x   - (lmSelect)
##   ... - forwarded to 'variable.names.lmSelect'
##
## Rval:  (formula)
##
formula.lmSelect <- function (x, ...) {
    ## lm
    if (is.null(x$.lm)) {
        stop ("'lmSelect' object does not have a 'lm' component")
    }

    ## full model
    f <- formula(x$.lm)
    ## variable names (complementary submodel)
    x.names <- variable.names(x, ..., .cmpl = TRUE)
    ## update formula
    f <- update(f, paste(c(".~.", x.names), collapse = "-"))

    ## done
    f
}


## extract model frame
##
## Args:
##   formula       - (lmSelect)
##   ...           - forwarded to 'formula.lmSelect'
##
## Rval:  (data.frame)
##
model.frame.lmSelect <- function (formula, ...) {
    ## lm
    if (is.null(formula$.lm)) {
        stop ("'lmSelect' object does not have a 'lm' component")
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
##   object        - (lmSelect)
##   best          - (integer)
##   ...           - ignored
##
## Rval:  (matrix)
##
model.matrix.lmSelect <- function (object, best = 1, ...) {
    ## lm
    if (is.null(object$.lm)) {
        stop ("'lmSelect' object does not have a 'lm' component")
    }

    ## full model
    x <- model.matrix(object$.lm)
    ## which
    wi <- object$which[, best]
    ## select
    x <- x[, wi]

    ## done
    x
}


## refit
##
## Args:
##   object - (lmSelect)
##   ...    - forwarded to 'formula.lmSelect'
##
## Rval:  (lm)
##
refit.lmSelect <- function (object, ...) {
    ## lm
    if (is.null(object$.lm)) {
        stop ("'lmSelect' object does not have a 'lm' component")
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
##   object - (lmSelect)
##   ...    - forwarded to 'refit.lmSelect'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   This method refits the 'lmSelect' object and calls 'coef' on the
##   obtained 'lm' object.  In some circumstances it might be more
##   efficient to call 'refit' directly and execute 'coef' on the
##   obtained 'lm' object.
##
coef.lmSelect <- function (object, ...) {
    coef(refit(object, ...))
}


## extract variance-covariance matrix
##
## Args:
##   object - (lmSelect)
##   ...    - forwarded to 'refit.lmSelect'
##
## Rval:  (matrix)
##
## Note: 'refit'
##   See 'coef.lmSelect'.
##
vcov.lmSelect <- function (object, ...) {
    vcov(refit(object, ...))
}


## extract fitted values
##
## Args:
##   object - (lmSelect)
##   ...    - forwarded to 'refit.lmSelect'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.lmSelect'.
##
fitted.lmSelect <- function (object, ...) {
  fitted(refit(object, ...))
}


## extract residuals
##
## Args:
##   object - (lmSelect)
##   ...    - forwarded to 'refit.lmSelect'
##
## Rval:  (numeric[])
##
## Note: 'refit'
##   See 'coef.lmSelect'.
##
residuals.lmSelect <- function (object, ...) {
  residuals(refit(object, ...))
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
    d <- object$rss[best]

    ## done
    d
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
        S <- sum(log(w))
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
##   k      - (numeric|"AIC"|"BIC")
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
