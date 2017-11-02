##
## File:  lmSelect.R
##



###########################
##  WORKHORSE FUNCTIONS  ##
###########################


## workhorse function
##
## Arguments:
##   x       - (double[,]) model matrix
##   y       - (double[]) response variable
##   weights - (double[]) 
##   offset  - (double[])
##   include - (integer[]) regressors to force in
##   exclude - (integer[]) regressors to force out
##   penalty - (character|numeric|function) penalty per parameter
##   nbest   - (integer) number of best subsets
##   pradius - (integer) preordering radius
##   ...     - ignored
##
## Result: (list)
##   NOBS         - (integer)
##   nobs         - (integer)
##   nvar         - (integer)
##   weights      - (double[])
##   intercept    - (logical)
##   include      - (integer[])
##   exclude      - (integer[])
##   sizes        - (integer[])
##   tolerance    - (double)
##   nbest        - (integer)
##   rank         - (integer[nbest])
##   df_residual  - (integer[nbest])
##   penalty      - (double[nbest])
##   which        - (logical[nvar,nbest])
##   .interrupted - (logical)
##   .nodes       - (integer)
##
lmSelect_fit <- function (x, y, weights = NULL, offset = NULL,
                          include = NULL, exclude = NULL,
                          penalty = "BIC", tolerance = 0,
                          nbest = 1, ..., pradius = NULL) {
    ## dims
    NOBS <- NROW(x)

    ## intercept
    has_intercept <- all(x[, 1] == 1)

    ## offset
    if (is.null(o <- offset))  o <- rep(0, NOBS)

    y <- y - o

    ## weights
    if (is.null(w <- weights))  w <- rep(1, NOBS)

    ok <- w != 0
    w <- sqrt(w[ok])

    x <- w * x[ok, , drop = FALSE]
    y <- w * y[ok,   drop = FALSE]

    ## dims
    nobs <- NROW(x)
    nvar <- NCOL(x)

    ## variables names
    x_names <- colnames(x)

    ## include
    if (is.null(include)) {
        include <- logical(nvar)
    } else {
        if (is.numeric(include)) {
            ## numeric --> numeric
            include <- match(include, seq_len(nvar))
        } else if (is.character(include)) {
            ## character --> numeric
            include <- match(include, x_names)
        } else if (is.logical(include)) {
            ## logical --> numeric
            include <- which(rep(include, length.out = nvar))
        }

        if (any(is.na(include))) {
            warning ("non-existing columns selected in 'include'")
        }

        ## canonicalize
        include <- is.element(seq_len(nvar), include)
    }

    names(include) <- x_names

    ## exclude
    if (is.null(exclude)) {
        exclude <- logical(nvar)
    } else {
        if (is.numeric(exclude)) {
            ## numeric --> numeric
            exclude <- match(exclude, seq_len(nvar))
        } else if (is.character(exclude)) {
            ## character --> numeric
            exclude <- match(exclude, x_names)
        } else if (is.logical(exclude)) {
            ## logical --> numeric
            exclude <- which(rep(exclude, length.out = nvar))
        }

        if (any(is.na(include))) {
            warning ("invalid columns selected in 'exclude'")
        }

        ## canonicalize
        exclude <- is.element(seq_len(nvar), exclude)
    }

    names(exclude) <- x_names

    ## include/exclude non-overlapping
    if (any(include & exclude)) {
        warning ("'include' and 'exclude' overlap: fixing 'exclude'")

        exclude <- exclude & !include
    }

    ## intercept
    if (has_intercept) {
        if (include[1]) {
            ## OK: included
        } else if (exclude[1]) {
            ## OK: excluded
        } else {
            ## default: include
            include[1] <- TRUE
        }
    }

    ## which
    which_include <- which(include)
    which_exclude <- which(exclude)

    which_maybe <- setdiff(seq_len(nvar), which_include)
    which_maybe <- setdiff(which_maybe, which_exclude)

    which_model <- c(which_include, which_maybe)

    ## size
    size <- length(which_model)
    mark <- length(which_include)
    nmin <- mark + 1
    nmax <- size

    ## preordering radius
    ## TODO: validate
    if (is.null(pradius)) {
        pradius <- floor((size - mark) / 3)
    }

    ## penalty
    pen <- get("penalty", envir = parent.env())
    pen <- do.call(pen, list(spec = penalty))

    ## call C
    z <- .Call("lmSelect", list(...)$.algo,
               cbind(x[, which_model], y), mark,
               arg_penalty(pen, nobs), tolerance + 1, nbest,
               pradius)

    ## result
    best_names <- format_ordinal(seq_len(nbest))

    ans <- list()
    ans$NOBS <- NOBS
    ans$nobs <- nobs
    ans$nvar <- nvar
    ans$weights <- weights
    ans$intercept <- has_intercept
    ans$include <- include
    ans$exclude <- exclude
    ans$sizes <- seq.int(nmin, nmax)
    ans$tolerance <- tolerance
    ans$nbest <- nbest

    ans$rank <- structure(z$rank, names = best_names)
    ans$df_residual <- ans$nobs - ans$rank

    ans$penalty <- structure(z$penalty, names = best_names, def = pen)

    ans$which <- array(NA, dim = c(nvar, nbest),
                       dimnames = list(x_names, best_names))
    ans$which[which_model, ] <- z$which

    ans$.interrupted <- z$.interrupted
    ans$.nodes <- z$.nodes

    ## done
    ans
}



#########################
##  GENERATOR METHODS  ##
#########################


## coercion method
##
## Arguments:
##   formula - (lmSelect)
##   ...     - ignored
##
## Result: (lmSelect)
##
lmSelect.lmSubsets <- function (formula, penalty = "BIC", ...) {
    x <- formula;  formula <- NULL

    ## tolerance
    if (!all(x$tolerance[x$sizes] == 0)) {
        stop ("non-zero tolerance")
    }

    ## call
    x$call[1] <- call("lmSelect")
    x$call$nmin <- x$call$nmax <- NULL
    x$call$penalty <- penalty

    ## penalty
    pen <- get("penalty", envir = parent.env())
    pen <- do.call(pen, list(spec = penalty))

    penalty <- eval_penalty(pen, x)

    ## order
    pi <- order(penalty)
    pi <- array(c((pi - 1) %%  x$nbest + 1,
                  (pi - 1) %/% x$nbest + 1),
                dim = c(length(pi), 2))
    nok <- is.na(penalty[pi])
    pi <- pi[!nok            , , drop = FALSE]
    pi <- pi[seq_len(x$nbest), , drop = FALSE]

    x_names <- dimnames(x$which)[[1]]
    best_names <- format_ordinal(seq_len(x$nbest))

    ## rank
    x$rank <- array(x$rank[pi], dim = x$nbest,
                    dimnames = list(best_names))

    ## penalty
    x$penalty <- array(penalty[pi], dim = x$nbest,
                       dimnames = list(best_names))
    attr(x$penalty, "def") <- pen

    ## which table
    which <- cbind(rep(seq_len(x$nvar), times = nrow(pi)),
                   rep(pi[, 1]        , each  = x$nvar),
                   rep(pi[, 2]        , each  = x$nvar))
    x$which <- array(x$which[which], dim = c(x$nvar, x$nbest),
                     dimnames = list(x_names, best_names))

    ## class
    class(x) <- "lmSelect"

    ## done
    x
}


## matrix interface
##
## Arguments:
##   formula  - (double[,])
##   y        - (double[])
##   interept - (logical)
##   ...      - forwarded
##
## Result: (lmSelect) see also 'lmSelect.default'
##
lmSelect.matrix <- function (formula, y, intercept = TRUE, ...) {
    x <- formula;  formula <- NULL

    ## names
    x_names <- paste0("x[, ", seq_len(ncol(x)), "]")
    if (all(x[, 1] == 1)) {
        x_names <- x_names[, -1]
        intercept <- TRUE
    }

    ## formula
    f <- stats_formula(x_names, "y", intercept)

    ## forward call
    cl <- match.call()
    cl[[1]] <- quote(lmSelect)
    cl$data <- bquote(list(x = .(x), y = .(y)), list(x = cl$formula, y = cl$y))
    cl$formula <- f
    cl$y <- cl$intercept <- NULL

    eval(cl, parent.frame())
}


## standard formula interface
##
## Arguments:
##   formula   - (formula) an object that can be coerced to class 'formula'
##   data      - (data.frame)
##   subset    - (integer[]) subset of observations
##   weights   - (double[])
##   na.action - (function)
##   model     - (logical)
##   x         - (logical)
##   y         - (logical)
##   contrasts - (list)
##   offset    - (double[])
##   ...       - forwarded to 'lmSelect_fit'
##
## Result: (lmSelect)  see 'lmSelect_fit'
##   call      - (call)
##   na.action - ()
##   offset    - (double[])
##   contrasts - (list)
##   xlevels   - (list)
##   terms     - (terms, formula)
##   model     - (data.frame) model frame
##   x         - (double[,])
##   y         - (double[])
##
lmSelect.default <- function (formula, data, subset, weights, na.action,
                              model = TRUE, x = FALSE, y = FALSE,
                              contrasts = NULL, offset, ...) {
    ## construct 'lmSelect' object
    ret_m <- model;  model <- NULL
    ret_x <- x;  x <- NULL
    ret_y <- y;  y <- NULL

    ## keep call (of generic)
    cl <- match.call()
    cl[[1]] <- as.name("lmSelect")

    ## model frame
    mf <- match.call(expand.dots = FALSE)
    m <- c("formula", "data", "subset",
           "weights", "na.action", "offset")
    m <- match(m, names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(model.frame)
    mf$drop.unused.levels <- TRUE
    mf <- eval(mf, parent.frame())

    ## model terms
    mt <- attr(mf, "terms")

    ## model matrix and response
    x <- model.matrix(mt, mf, contrasts)
    y <- model_response(mf, "numeric")

    ## weights
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) {
        stop ("'weights' must be a numeric vector")
    }

    ## offset
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset) && (length(offset) != NROW(y))) {
        stop ("length of 'offset' must equal number of observations")
    }

    ## fit subsets
    ans <- lmSelect_fit(x, y, w, offset = offset, ...)

    ## result
    class(ans) <- "lmSelect"

    ans$call <- cl
    ans$na.action <- attr(mf, "na.action")
    ans$offset <- offset
    ans$contrasts <- attr(x, "contrasts")
    ans$xlevels <- .getXlevels(mt, mf)
    ans$terms <- mt

    if (ret_m)  ans$model <- mf
    if (ret_x)  ans$x <- x
    if (ret_y)  ans$y <- y

    ## done
    ans
}



########################
##  STANDARD METHODS  ##
########################


## format method for 'lmSelect' objects
##
## Arguments:
##   x   - (lmSelect)
##   ... - forwarded
##
## Result: (list)
##
format.lmSelect <- function (x, ...) {
    ans <- list()

    ## call
    ans$call <- paste(deparse(x$call), sep = "\n", collapse = "\n")

    ## penalty
    ans$penalty <- ifelse(is.na(x$penalty), "",
                          format_numeric(x$penalty, ...))

    ## which
    ans$which <- apply(x$which[, , drop = FALSE], 1, format_which)
    ans$which <- matrix(ans$which, dimnames = list(names(ans$which), "best"))
    rownames(ans$which)[x$include] <- paste0("+", rownames(ans$which)[x$include])
    rownames(ans$which)[x$exclude] <- paste0("-", rownames(ans$which)[x$exclude])

    ## done
    ans
}


## print method for 'lmSelect' objects
##
## Arguments:
##   x   - (lmSelect)
##   ... - forwarded
##
## Result: (lmSelect) invisible
##
print.lmSelect <- function (x, ...) {
    catln <- function (...)  base::cat(..., "\n", sep = "")
    print <- function (..., skip = 0, indent = 0) {
        output <- capture.output(base::print(...))
        if (skip > 0) output <- output[-seq_len(skip)]
        indent <- paste0(rep(" ", indent), collapse = "")
        cat(paste0(indent, output, "\n"), sep = "")
    }

    ## format
    fmt <- format(x, ...)

    catln()

    ## call
    catln("Call:")
    catln("  ", fmt$call)

    catln()

    ## fit
    catln("Model fit:")
    catln("  best")
    catln("    ", format_penalty(attr(x$penalty, "def")))
    print(fmt$penalty, quote = FALSE, indent = 2)

    catln()

    ## which
    catln("Which:")
    print(fmt$which, quote = FALSE, indent = 2)

    catln()

    ## done
    invisible(x)
}


## plot method for 'lmSelect' objects
##
## Arguments:
##   x    - (lmSelect)
##   y    - ignored
##   ...  - forwarded
##   axes - (logical) draw axes
##   ann  - (logical) annotate
##
## Result: (lmSelect) invisible
##
plot.lmSelect <- function (x, ..., axes = TRUE, ann = par("ann")) {
    object <- x;  x <- NULL

    ## new plot
    plot.new()

    ## pars
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))

    ## margins
    par(mar = c(5, 4, 4, 4) + 0.1)

    ## penalty
    xlim <- c(1, object$nbest)
    ylim <- range(object$penalty, na.rm = TRUE)
    plot.window(xlim, ylim, ...)

    lines(seq_len(object$nbest), object$penalty, type = "o", lty = 1,
          pch = 16, col = "red", bg = "white", ...)

    ## axes
    if (axes) {
        axis(1)
        axis(2)
        box()
    }

    ## annotations
    if (ann) {
        title(main = "Best subsets", xlab = "best", ylab = "Penalty")

        legend <- format_penalty(attr(object$penalty, "def"))
        legend("topleft", legend = paste0("Penalty: ", legend),
               lty = 1, pch = 16, col = "red", pt.bg = "white",
               bty = "n")
    }

    ## done
    invisible(x)
}



#########################
##  EXTRACTOR METHODS  ##
#########################


## is NA
##
## Arguments:
##   x    - (lmSelect)
##   best - (integer[])
##   ...  - ignored  
##
## Result: (logical[])
##
is_NA.lmSelect <- function (x, best = 1, ...) {
    if (is.null(best)) {
        best <- seq_len(x$nbest)
    }

    is.na(x$rank[best])
}


## extract variable names
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##
## Result: (character[]) variable names
##
variable.names.lmSelect <- function (object, best = 1, ...) {
    ## full model
    x_names <- dimnames(object$which)[[1]]

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, best)) {
        return (NA)
    }

    ## submodel
    which <- object$which[, best]
    x_names[which]
}


## extract formula
##
## Arguments:
##   x    - (lmSelect)
##   best - (integer)
##   ...  - ignored
##
## Result: (formula)
##
formula.lmSelect <- function (x, best, ...) {
    ## 'best' processing
    if (missing(best)) {
        f <- formula(x$terms)

        return (f)
    }

    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(x, best)) {
        return (NA)
    }

    ## names
    y_name <- all.vars(x$terms)[1]
    x_names <- variable.names(x, best)

    if (x$intercept) {
        x_names <- x_names[-1]
    }

    ## build formula
    stats_formula(x_names, y_name, x$intercept,
                  env = environment(x$terms))
}


## extract model frame
##
## Arguments:
##   object - (lmSubsets)
##   best   - (integer)
##   ...    - ignored
##
## Result: (character[]) variable names
##
model.frame.lmSelect <- function(formula, best, ...) {
    ## further arguments
    args <- list(...)
    m <- c("data", "na.action", "subset")
    m <- match(m, names(args), 0L)
    args <- args[m]
    
    ## 'best' processing
    if (missing(best)) {
        if ((length(args) == 0) && !is.null(mf <- formula[["model"]]))  return (mf)
    } else if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (!missing(best) && is_NA(formula, best)) {
        return (NA)
    }

    ## extract formula
    if (!missing(best)) {
        f <- formula(formula, best = best)
    } else {
        f <- terms(formula)
    }

    ## forward call to 'model.frame'
    cl <- formula$call
    m <- c("formula", "data", "subset",
           "weights", "na.action", "offset")
    m <- match(m, names(cl), 0L)
    cl <- cl[c(1L, m)]
    cl[[1L]] <- quote(model.frame)
    cl$formula <- f
    cl$drop.unused.levels <- TRUE
    cl$xlev <- formula$xlevels

    env <- environment(formula$terms)
    if (is.null(env)) {
        env <- parent.frame()
    }

    eval(cl, env)
}


## extract model matrix
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - forwarded
##
## Result: (double[,])
##
model.matrix.lmSelect <- function (object, best, ...) {
    ## full model
    x <- object[["x"]]
    if (is.null(x)) {
        data <- model.frame(object, xlev = object$xlevels, ...)
        x <- NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
    }

    ## 'best' processing
    if (missing(best)) {
        return (x)
    }

    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, best)) {
        return (NA)
    }

    ## submodel
    which <- object$which[, best]
    x[, which]
}


## extract model response
##
## Arguments:
##   x   - (lmSubsets)
##   ... - ignored
##
## Result: (formula)
##
model_response.lmSelect <- function (data, ...) {
    ## extract response
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
## Arguments:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - forwarded
##
## Result: (lm)
##
refit.lmSelect <- function (object, best = 1, ...) {
    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, best)) {
        return (NA)
    }

    ## extract formula
    f <- formula(object, best = best)

    ## forward call to 'lm'
    cl <- object$call
    m <- c("formula", "data", "subset",
           "weights", "na.action", "offset")
    m <- match(m, names(cl), 0L)
    cl <- cl[c(1L, m)]
    cl[[1L]] <- quote(lm)
    cl$formula <- f

    dots <- list(...)
    cl[names(dots)] <- dots

    env <- environment(object$terms)
    if (is.null(env)) {
        env <- parent.frame()
    }

    eval(cl, env)
}


## extract deviance
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer[])
##   ...    - ignored
##
## Result: (double[])
##
deviance.lmSelect <- function (object, best = 1, ...) {
    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(object$nbest)
    }

    ## weights
    if (is.null(w <- object$weights)) {
        w <- rep(1, object$NOBS)
    } else {
        w <- sqrt(w)
    }

    ## extract RSS
    structure(
        sapply(best, function (i) {
            r <- residuals(object, best = i)
            sum((w * r)^2)
        }),
        names = format_ordinal(best)
    )
}


## extract residual standard deviation
##
## Arguments:
##   object - (lmSubsets)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double)
##
sigma.lmSelect <- function (object, best = 1, ...) {
    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, best)) {
        return (NA)
    }

    ## compute
    rss <- deviance(object, best = best)
    rdf <- object$nobs - object$rank[best]
    sigma2 <- stats_sigma2(rss, rdf)
    sqrt(sigma2)
}


## extract log-likelihood
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer[])
##   ...    - ignored
##   drop   - (logical)
##
## Result: (double[])
##
log_lik.lmSelect <- function (object, best = 1, ...) {
    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(object$nbest)
    }

    ## extract rss
    rss <- deviance(object, best = best)

    ## compute log-likelihoods
    ans <- stats_log_lik(object$nobs, object$weights, rss)

    ## done
    ans
}


## extract log-likelihood
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##
## Result: (logLik)
##
logLik.lmSelect <- function (object, best = 1, ...) {
    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, best)) {
        return (NA)
    }

    ## compute log-likelihoods
    ans <- log_lik(object, best = best)

    ## decorate
    ans <- structure(ans, df = object$rank[best] + 1,
                     nobs = object$nobs)
    class(ans) <- "logLik"

    ## done
    ans
}


## compute AIC
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer[])
##   ...    - ignored
##   k      - (double) penalty
##
## Result: (double[])
##
aic.lmSelect <- function (object, best = 1, ..., k = 2) {
    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(object$nbest)
    }

    ## extract log-likelihoods
    ll <- log_lik(object, best = best)

    ## extract DF
    df <- object$rank[best] + 1

    ## compute AICs
    ans <- stats_aic(ll, k, df)

    ## done
    ans
}


## compute AIC
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer[])
##   ...    - ignored
##   k      - (double) penalty
##   drop   - (logical)
##
## Result: (double|data.frame)
##
AIC.lmSelect <- function (object, best = 1, ..., k = 2) {
    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(object$nbest)
    }

    ## compute AIC
    ans <- aic(object, best = best)

    ## decorate
    if (length(ans) > 1) {
        ans <- data.frame(df = object$rank[best] + 1,
                          AIC = ans)
    }

    ## done
    ans
}


## compute BIC
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer[])
##   ...    - ignored
##
## Result: (double[])
##
bic.lmSelect <- function (object, best = 1, ...) {
    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(object$nbest)
    }

    ## extract log-likelihoods
    ll <- log_lik(object, best = best)

    ## extract DF
    df <- object$rank[best] + 1

    ## compute AICs
    ans <- stats_bic(ll, object$nobs, df)

    ## done
    ans
}


## compute BIC
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer[])
##   ...    - ignored
##
## Result: (double|data.frame)
##
BIC.lmSelect <- function (object, best = 1, ...) {
    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(object$nbest)
    }

    ## compute BIC
    ans <- bic(object, best = best)

    ## decorate
    if (length(ans) > 1) {
        ans <- data.frame(df = object$rank[best] + 1,
                          BIC = ans)
    }

    ## done
    ans
}


## extract ceofficients
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double[])
##
coef.lmSelect <- function (object, best = 1, ...) {
    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, best)) {
        return (NA)
    }

    ## extract submodel
    x <- model.matrix(object, best = best)
    y <- model_response(object)

    ## solve
    ans <- stats_coef(x, y, object$offset, object$weights)

    ## names
    x_names <- variable.names(object, best = best)
    names(ans) <- x_names

    ## done
    ans
}


## extract variance-covariance matrix
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double[,])
##
vcov.lmSelect <- function (object, best = 1, ...) {
    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, best)) {
        return (NA)
    }

    ## extract submodel
    x <- model.matrix(object, best = best)
    y <- model_response(object)

    ## solve
    ans <- stats_vcov(x, y, object$offset, object$weights)

    ## names
    x_names <- variable.names(object, best = best)
    dimnames(ans) <- list(x_names, x_names)

    ## done
    ans
}


## extract fitted values
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double[])
##
fitted.lmSelect <- function (object, best = 1, ...) {
    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, best)) {
        return (NA)
    }

    ## extract submodel
    x <- model.matrix(object, best = best)
    y <- model_response(object)

    ## solve
    ans <- stats_fitted(x, y, object$offset, object$weights)

    ## NA action
    ans <- napredict(object$na.action, ans)

    ## done
    ans
}


## extract residuals
##
## Arguments:
##   object - (lmSelect)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double[])
##
residuals.lmSelect <- function (object, best = 1, ...) {
    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, best)) {
        return (NA)
    }

    ## extract submodel
    x <- model.matrix(object, best = best)
    y <- model_response(object)

    ## solve
    ans <- stats_residuals(x, y, object$offset, object$weights)

    ## NA action
    ans <- naresid(object$na.action, ans)

    ## done
    ans
}
