##
## File:  lmSubsets.R
##



###########################
##  WORKHORSE FUNCTIONS  ##
###########################


## workhorse function
##
## Arguments:
##   x         - (double[,]) model matrix
##   y         - (double[]) model response
##   weights   - (double[]) 
##   offset    - (double[])
##   include   - (integer[]) regressors to force in
##   exclude   - (integer[]) regressors to force out
##   nmin      - (integer)
##   nmax      - (integer)
##   tolerance - (double[])
##   nbest     - (integer) number of best subsets
##   ...       - ignored
##   pradius   - (integer) preordering radius
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
##   tolerance    - (double[])
##   nbest        - (integer)
##   rank         - (integer[nbest,nvar])
##   df_residual  - (integer[nbest,nvar])
##   rss          - (double[nbest,nvar])
##   which        - (logical[nvar,nbest,nvar])
##   .interrupted - (logical)
##   .nodes       - (integer)
##
lmSubsets_fit <- function (x, y, weights = NULL, offset = NULL,
                           include = NULL, exclude = NULL,
                           nmin = NULL, nmax = NULL, tolerance = 0,
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
    nmin <- max(nmin, mark + 1)
    nmax <- min(nmax, size)

    ## preordering radius
    ## TODO: validate
    if (is.null(pradius)) {
        pradius <- floor((size - mark) / 3)
    }

    ## tolerance
    ## TODO: validate
    tolerance <- rep(tolerance, length.out = size)
    tolerance[-seq.int(nmin, nmax)] <- +Inf

    ## call C
    z <- .Call("lmSubsets", list(...)$.algo,
               cbind(x[, which_model], y), mark, tolerance + 1,
               nbest, pradius)

    ## result
    best_names <- format_ordinal(seq_len(nbest))
    size_names <- seq_len(nvar)

    ans <- list()
    ans$NOBS <- NOBS
    ans$nobs <- nobs
    ans$nvar <- nvar
    ans$weights <- weights
    ans$intercept <- has_intercept
    ans$include <- include
    ans$exclude <- exclude
    ans$sizes <- seq.int(nmin, nmax)
    ans$tolerance <- array(tolerance, dim = nvar,
                           dimnames = list(size_names))
    ans$nbest <- nbest

    ans$rank <- array(NA, dim = c(nbest, nvar),
                      dimnames = list(best_names, size_names))
    ans$rank[, seq_len(size)] <- z$rank
    ans$df_residual <- ans$nobs - ans$rank

    ans$rss <- array(NA, dim = c(nbest, nvar),
                     dimnames = list(best_names, size_names))
    ans$rss[, seq_len(size)] <- z$rss

    ans$which <- array(NA, dim = c(nvar, nbest, nvar),
                       dimnames = list(x_names, best_names, size_names))
    ans$which[which_model, , seq_len(size)] <- z$which

    ans$.interrupted <- z$.interrupted
    ans$.nodes <- z$.nodes

    ## done
    ans
}



#########################
##  GENERATOR METHODS  ##
#########################


## matrix interface
##
## Arguments:
##   formula  - (double[,])
##   y        - (double[])
##   interept - (logical)
##   ...      - forwarded
##
## Result: (lmSubsets) see also 'lmSubsets.default'
##
lmSubsets.matrix <- function (formula, y, intercept = TRUE, ...) {
    x <- formula;  formula <- NULL

    ## names
    x_names <- paste0("x[, ", seq_len(ncol(x)), "]")
    if (all(x[, 1] == 1)) {
        x_names <- x_names[-1]
        intercept <- TRUE
    }

    ## formula
    f <- stats_formula(x_names, "y", intercept)

    ## forward call
    cl <- match.call()
    cl[[1]] <- quote(lmSubsets)
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
##   subset    - (vector) subset of observations
##   weights   - (double[])
##   na.action - (function)
##   model     - (logical)
##   x         - (logical)
##   y         - (logical)
##   contrasts - (list)
##   offset    - (double[])
##   ...       - forwarded to 'lmSubsets_fit'
##
## Result: (lmSubsets) see also 'lmSubsets_fit'
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
lmSubsets.default <- function (formula, data, subset, weights, na.action,
                               model = TRUE, x = FALSE, y = FALSE,
                               contrasts = NULL, offset, ...) {
    ret_m <- model;  model <- NULL
    ret_x <- x;  x <- NULL
    ret_y <- y;  y <- NULL

    ## keep call (of generic)
    cl <- match.call()
    cl[[1]] <- as.name("lmSubsets")

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
    ans <- lmSubsets_fit(x, y, weights = w, offset = offset, ...)

    ## result
    class(ans) <- "lmSubsets"

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


## format method for 'lmSubsets' objects
##
## Arguments:
##   x   - (lmSubsets)
##   ... - forwarded
##
## Result: (list)
##
format.lmSubsets <- function (x, ...) {
    ans <- list()

    ## call
    ans$call <- paste(deparse(x$call), sep = "\n", collapse = "\n")

    ## rss
    ans$rss <- x$rss[, x$sizes, drop = FALSE]
    ans$rss <- ifelse(is.na(ans$rss), "", format_numeric(ans$rss, ...))
    colnames(ans$rss) <- paste0(colnames(ans$rss), " (",
                                format_numeric(x$tolerance[x$sizes], ...),
                                ")")

    ## which
    ans$which <- apply(x$which[, , , drop = FALSE], c(1, 2), format_which)
    rownames(ans$which)[x$include] <- paste0("+", rownames(ans$which)[x$include])
    rownames(ans$which)[x$exclude] <- paste0("-", rownames(ans$which)[x$exclude])

    ## done
    ans
}


## print method for 'lmSubsets' objects
##
## Arguments:
##   x       - (lmSubsets)
##   ...     - forwarded
##
## Result: (lmSubsets) invisible
##
print.lmSubsets <- function (x, ...) {
    catln <- function (...)  base::cat(..., "\n", sep = "")
    print <- function (..., skip = 0, indent = 0) {
        output <- capture.output(base::print(...))
        if (skip > 0) output <- output[-seq_len(skip)]
        indent <- paste0(rep(" ", indent), collapse = "")
        cat(paste0(indent, output, "\n"), sep = "")
    }

    ## format
    fmt <- format(x, digits = 3, ...)

    catln()

    ## call
    catln("Call:")
    catln("  ", fmt$call)

    catln()

    ## fit
    catln("Model fit:")
    catln("  best x size (tolerance)")
    catln("    deviance")
    print(fmt$rss, quote = FALSE, indent = 2)

    catln()

    ## which
    catln("Which:")
    catln("  variable x best")
    catln("    sizes")
    print(fmt$which, quote = FALSE, indent = 2)

    catln()

    ## done
    invisible(x)
}


## plot method for 'lmSubsets' objects
##
## Arguments:
##   x       - (lmSubsets)
##   penalty - (character|numeric|function)
##   ...     - forwarded
##   axes    - (logical) draw axes
##   ann     - (logical) annotate
##
## Result: (lmSubsets) invisible
##
plot.lmSubsets <- function (x, penalty = "BIC", ...,
                            axes = TRUE, ann = par("ann")) {
    object <- x;  x <- NULL

    ## new plot
    plot.new()

    ## pars
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))

    ## margins
    par(mar = c(5, 4, 4, 4) + 0.1)

    ## RSS
    xlim <- range(object$sizes)
    ylim <- range(object$rss, na.rm = TRUE)
    plot.window(xlim, ylim, ...)

    for (j in object$sizes) {
        x <- rep_len(j, object$nbest)
        y <- rev(object$rss[seq_len(object$nbest), j])

        lines(x, y, type = "o", lty = 3,
              pch = c(rep_len(21, object$nbest - 1), NA),
              col = "black", bg = "white", ...)
    }

    lines(object$sizes, object$rss[1, object$sizes], type = "o",
          lty = 1, pch = 16, col = "black", bg = "white", ...)

    ## axes
    if (axes) {
        axis(1)
        axis(2)
        box()
    }

    ## annotations
    if (ann) {
        title(main = "All subsets",
              sub = paste("nbest = ", object$nbest),
              xlab = "Number of regressors", ylab = "Deviance")

        legend <- "Deviance: RSS"
    }

    ## penalty
    if (!is.null(penalty)) {
        pen <- get("penalty", envir = parent.env())
        pen <- do.call(pen, list(spec = penalty))

        val <- eval_penalty(pen, object)

        ylim <- range(val, na.rm = TRUE)
        plot.window(xlim, ylim, ...)

        for (j in object$sizes) {
            x <- rep_len(j, object$nbest)
            y <- rev(val[seq_len(object$nbest), j])

            lines(x, y, type = "o", lty = 3,
                  pch = c(rep_len(21, object$nbest - 1), NA),
                  col = "red", bg = "white", ...)
        }

        lines(object$sizes, val[1, object$sizes], type = "o", lty = 1,
              pch = 16, col = "red", bg = "white", ...)

        ## axes
        if (axes) {
            axis(4)
        }

        ## annotations
        if (ann) {
            mtext("Penalty", side = 4, line = 3)

            legend <- c(legend, paste0("Penalty: ", format_penalty(pen)))
        }
    }

    ## legend
    if (ann) {
        legend("topright", legend = legend, lty = 1, pch = 16,
               col = c("black", "red"), pt.bg = "white", bty = "n")
    }

    ## done
    invisible(object)
}



#########################
##  EXTRACTOR METHODS  ##
#########################


## is NA
##
## Arguments:
##   x    - (lmSubsets)
##   size - (integer[])
##   best - (integer[])
##   ...  - ignored
##   drop - (logical)
##
## Result: (logical[,])
##
is_NA.lmSubsets <- function (x, size, best = 1, ..., drop = TRUE) {
    ## 'size' processing
    if (missing(size)) {
        size <- x$sizes
    } else if (is.null(size)) {
        size <- seq_len(x$nvar)
    }

    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(x$nbest)
    }

    ans <- is.na(x$rank[best, size, drop = FALSE])
    if (drop)  ans <- drop(ans)

    ## done
    ans
}


## extract variable names
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Result: (character[])
##
variable.names.lmSubsets <- function (object, size, best = 1, ...) {
    ## full model
    x_names <- dimnames(object$which)[[1]]

    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, size, best)) {
        return (NA)
    }

    ## submodel
    which <- object$which[, best, size]
    x_names[which]
}


## extract formula
##
## Arguments:
##   x    - (lmSubsets)
##   size - (integer)
##   best - (integer)
##   ...  - ignored
##
## Result: (formula)
##
formula.lmSubsets <- function (x, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        f <- formula(x$terms)

        return (f)
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(x, size, best)) {
        return (NA)
    }

    ## names
    y_name <- all.vars(x$terms)[1]
    x_names <- variable.names(x, size = size, best = best)

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
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Result: (model.frame)
##
model.frame.lmSubsets <- function (formula, size, best = 1, ...) {
    ## further arguments
    args <- list(...)
    m <- c("data", "na.action", "subset")
    m <- match(m, names(args), 0L)
    args <- args[m]
    
    ## 'size' processing
    if (missing(size)) {
        if ((length(args) == 0) && !is.null(mf <- formula[["model"]]))  return (mf)
    } else if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (!missing(size) && is_NA(formula, size, best)) {
        return (NA)
    }

    ## extract formula
    if (!missing(size)) {
        f <- formula(formula, size = size, best = best)
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
    cl[names(args)] <- args

    env <- environment(formula$terms)
    if (is.null(env)) {
        env <- parent.frame()
    }

    eval(cl, env)
}


## extract model matrix
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - forwarded
##
## Result: (double[,])
##
model.matrix.lmSubsets <- function (object, size, best = 1, ...) {
    ## full model
    x <- object[["x"]]
    if (is.null(x)) {
        data <- model.frame(object, xlev = object$xlevels, ...)
        x <- NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
    }

    ## 'size' processing
    if (missing(size)) {
        return (x)
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, size, best)) {
        return (NA)
    }

    ## submodel
    which <- object$which[, best, size]
    x[, which]
}


## extract model response
##
## Arguments:
##   data - (lmSubsets)
##   ...  - ignored
##
## Result: (formula)
##
model_response.lmSubsets <- function (data, ...) {
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
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - forwarded
##
## Result: (lm)
##
refit.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing (size)) {
        stop ("missing argument: 'size'")
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, size, best)) {
        return (NA)
    }

    ## extract formula
    f <- formula(object, size = size, best = best)

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
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##   drop   - (logical)
##
## Result: (double[,])
##
deviance.lmSubsets <- function (object, size, best = 1, ..., drop = TRUE) {
    ## 'size' processing
    if (missing(size)) {
        size <- object$sizes
    } else if (is.null(size)) {
        size <- seq_len(object$nvar)
    }

    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(object$nbest)
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


## extract residual standard deviation
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double)
##
sigma.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, size, best)) {
        return (NA)
    }

    ## compute
    rss <- deviance(object, size = size, best = best)
    sigma2 <- stats_sigma2(rss, object$nobs - size)
    sqrt(sigma2)
}


## extract log-likelihood
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##   drop   - (logical)
##
## Result: (double[,])
##
log_lik.lmSubsets <- function (object, size, best = 1, ..., drop = TRUE) {
    ## 'size' processing
    if (missing(size)) {
        size <- object$sizes
    } else if (is.null(size)) {
        size <- seq_len(object$nvar)
    }

    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(object$nbest)
    }

    ## extract RSS
    rss <- object$rss[best, size, drop = FALSE]

    ## compute log-likelihoods
    ans <- stats_log_lik(object$nobs, object$weights, rss)
    ans <- t(ans)

    ## drop
    if (drop) {
        ans <- drop(ans)
    }

    ## done
    ans
}


## extract log-likelihood
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Result: (logLik)
##
logLik.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing (size)) {
        stop ("missing argument: 'size'")
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, size, best)) {
        return (NA)
    }

    ## compute log-likelihoods
    ans <- log_lik(object, size = size, best = best)

    ## decorate
    ans <- structure(ans, df = size + 1,
                     nobs = object$ nobs)
    class(ans) <- "logLik"

    ## done
    ans
}


## compute AIC
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##   k      - (double)
##   drop   - (logical)
##
## Result: (double[,])
##
aic.lmSubsets <- function (object, size, best = 1, ..., k = 2, drop = TRUE) {
    ## 'size' processing
    if (missing(size)) {
        size <- object$sizes
    } else if (is.null(size)) {
        size <- seq_len(object$nvar)
    }

    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(object$nbest)
    }

    ## extract log-likelihood
    ll <- log_lik(object, size = size, best = best, drop = FALSE)

    ## extract DF
    df <- object$rank[best, size, drop = FALSE] + 1
    df <- t(df)

    ## compute AIC
    ans <- stats_aic(ll, k, df)

    ## drop
    if (drop) {
        ans <- drop(ans)
    }

    ## done
    ans
}


## compute AIC
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer)
##   ...    - ignored
##   k      - (double) penalty
##
## Result: (double|data.frame)
##
AIC.lmSubsets <- function (object, size, best = 1, ..., k = 2) {
    ## 'size' processing
    if (missing(size)) {
        size <- object$sizes
    } else if (is.null(size)) {
        size <- seq_len(object$nvar)
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## compute AIC
    ans <- aic(object, size = size, best = best, k = k)

    ## decorate
    if (length(ans) > 1) {
        ans <- data.frame(df = object$rank[best, size] + 1,
                          AIC = ans)
    }

    ## done
    ans
}


## compute BIC
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##   drop   - (logical)
##
## Result: (double[,])
##
bic.lmSubsets <- function (object, size, best = 1, ..., drop = TRUE) {
    ## 'size' processing
    if (missing(size)) {
        size <- object$sizes
    } else if (is.null(size)) {
        size <- seq_len(object$nvar)
    }

    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(object$nbest)
    }

    ## extract log-likelihood
    ll <- log_lik(object, size = size, best = best, drop = FALSE)

    ## extract DF
    df <- object$rank[best, size, drop = FALSE] + 1
    df <- t(df)

    ## compute AIC
    ans <- stats_bic(ll, object$nobs, df)

    ## drop
    if (drop) {
        ans <- drop(ans)
    }

    ## done
    ans
}


## compute BIC
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##
## Result: (double)
##
BIC.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        size <- object$sizes
    } else if (is.null(size)) {
        size <- seq_len(object$nvar)
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## compute BIC
    ans <- bic(object, size = size, best = best)

    ## decorate
    if (length(ans) > 1) {
        ans <- data.frame(df = object$rank[best, size] + 1,
                          BIC = ans)
    }

    ## done
    ans
}


## extract ceofficients
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double[])
##
coef.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, size, best)) {
        return (NA)
    }

    ## extract submodel
    x <- model.matrix(object, size = size, best = best)
    y <- model_response(object)

    ## solve
    ans <- stats_coef(x, y, object$offset, object$weights)

    ## names
    x_names <- variable.names(object, size = size, best = best)
    names(ans) <- x_names

    ## done
    ans
}


## extract variance-covariance matrix
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double[,])
##
vcov.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, size, best)) {
        return (NA)
    }

    ## extract submodel
    x <- model.matrix(object, size = size, best = best)
    y <- model_response(object)

    ## solve
    ans <- stats_vcov(x, y, object$offset, object$weights)

    ## names
    x_names <- variable.names(object, size = size, best = best)
    dimnames(ans) <- list(x_names, x_names)

    ## done
    ans
}


## extract fitted values
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double[])
##
fitted.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, size, best)) {
        return (NA)
    }

    ## extract submodel
    x <- model.matrix(object, size = size, best = best)
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
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double[])
##
residuals.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## available?
    if (is_NA(object, size, best)) {
        return (NA)
    }

    ## extract submodel
    x <- model.matrix(object, size = size, best = best)
    y <- model_response(object)

    ## solve
    ans <- stats_residuals(x, y, object$offset, object$weights)

    ## NA action
    ans <- naresid(object$na.action, ans)

    ## done
    ans
}
