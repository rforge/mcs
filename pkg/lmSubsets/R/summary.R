##
## File:  summary.R
##



#############
## SUMMARY ##
#############


## summary for 'lmSubsets' objects
##
## Args:
##   object  - (lmSubsets)
##   penalty - ("AIC"|"BIC"|numeric)
##   ...     - ignored
##
## Rval: (summary.lmSubsets)
##
summary.lmSubsets <- function (object, penalty = "AIC", ...) {

    ## aic
    N <- object$nobs
    if (is.null(w <- object$weights)) {
        S <- 0
    } else {
        S <- sum(log(w))
    }

    if (tolower(penalty) == "aic") {
        penalty <- structure("AIC", k = 2)
    } else if (tolower(penalty) == "bic") {
        penalty <- structure("BIC", k = log(N))
    } else if (is.numeric(penalty)) {
        penalty <- structure("AIC", k = penalty)
    } else {
        stop ("invalid 'penalty'")
    }

    ll <- 0.5 * (S - N * (log(2 * pi) + 1 - log(N) + log(object$rss)))
    aic <- -2 * ll + attr(penalty, "k") * object$df

    object$summary <- list(penalty = penalty, aic = aic)

    ## class
    class(object) <- c("summary.lmSubsets", "lmSubsets")

    ## done
    object
}


## summary for 'lmSelect' objects
##
## Args:
##   object  - (mcsSubset)
##   penalty - (("AIC"|"BIC"|numeric)[])
##   ...     - ignored
##
## Rval: (summary.mcsSubset)
##
summary.lmSelect <- function (object, penalty = "BIC", ...) {
    ## log lik
    N <- object$nobs
    if (is.null(w <- object$weights)) {
        S <- 0
    } else {
        S <- sum(log(w))
    }

    ll <- 0.5 * (S - N * (log(2 * pi) + 1 - log(N) + log(object$rss)))

    ## penalty
    kk <- NULL
    pp <- NULL
    for (p in penalty) {
        if (tolower(p) == "aic") {
            kk <- c(kk, 2    )
            pp <- c(pp, "AIC")
        } else if (tolower(p) == "bic") {
            kk <- c(kk, log(N))
            pp <- c(pp, "BIC" )
        } else if (is.numeric(p)) {
            kk <- c(kk, p    )
            pp <- c(pp, "AIC")
        } else {
            stop ("invalid 'penalty'")
        }
    }

    ## AIC
    aic <- matrix(rep(-2, length(kk)), ncol = 1) %*% ll        +
           matrix(kk                 , ncol = 1) %*% object$df

    ## summary
    object$summary <- list(penalty = structure(pp, k = kk),
                           aic     = aic                  )

    ## class
    class(object) <- c("summary.lmSelect", "lmSelect")

    ## done
    object
}



###########
## PRINT ##
###########


## print 'lmSubsets' summary
##
## Args:
##   x      - (summary.lmSubsets)
##   ...    - ignored
##
## Rval: (summary.lmSubsets) invisible
##
print.summary.lmSubsets <- function (x, ...)
{
    catln <- function (..., sep = "") base::cat(..., "\n", sep = sep)
    paste <- function (..., sep = "") base::paste(..., sep = sep)


    ## call
    catln()
    catln("Call:")
    catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))


    ## arguments
    catln()
    cat("Arguments:")
    val <- as.matrix(c(paste(x$summary$penalty, " (k = ", format(attr(x$summary$penalty, "k"), nsmall = 2), ")")),
                     ncol = 1)
    colnames(val) <- ""
    rownames(val) <- paste("  ", c("Value"), ":")
    print(val, quote = FALSE)


    ## fit
    catln()
    catln("Model fit (value):")
    catln("  best x size")
    aic <- x$summary$aic[, x$size, drop = FALSE]
    fit <- ifelse(is.na(aic), "", format(aic, nsmall = 2))
    rownames(fit) <- paste("  ", rownames(fit))
    print(fit, quote = FALSE)
    catln()


    ## done
    invisible(x)
}


## print 'lmSelect' summary
##
## Args:
##   x      - (summary.lmSelect)
##   digits - (integer)
##   ...    - ignored
##
## Rval: (summary.lmSelect) invisible
##
print.summary.lmSelect <- function (x, ...)
{
    catln <- function (..., sep = "") base::cat(..., "\n", sep = sep)
    paste <- function (..., sep = "") base::paste(..., sep = sep)


    ## call
    catln()
    catln("Call:")
    catln("  ", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))


    ## arguments
    catln()
    cat("Arguments:")
    val <- as.matrix(c(paste(x$summary$penalty, paste("(k = ", format(attr(x$summary$penalty, "k"), nsmall = 2), ")"), sep = " ", collapse = ", ")),
                     ncol = 1)
    colnames(val) <- ""
    rownames(val) <- paste("  ", c("Value"), ":")
    print(val, quote = FALSE)


    ## fit
    catln()
    catln("Model fit:")
    fit <- format(rbind(x$df, x$rss, x$summary$aic), nsmall = 2)
    rownames(fit) <- paste("  ", c("df", "Deviance", "Value", rep("", length(x$summary$penalty) - 1)))
    print(fit, quote = FALSE)
    catln()


    ## done
    invisible(x)
}



##########
## PLOT ##
##########


## plot 'lmSubsets' summary
##
## Args:
##   x      - (lmSubsets)
##   ...    - passed to 'plot.lmSubsets' and to 'matplot'.
##   legend - (character[])
##
## Rval: (summary.lmSubsets) invisible
##
plot.summary.lmSubsets <- function (x, ..., legend) {
    localPlot <- function (object, main, sub = NULL, xlab, ylab,
                           type, lty, pch, col, bg, ...) {
        if (missing(main)) main <- "All subsets (summary)"
        if (missing(xlab)) xlab <- "Number of regressors"
        if (missing(ylab)) ylab <- c("Deviance", "Value")

        type <- if (missing(type)) c("o", "o")         else rep(type, length = 2)
        lty  <- if (missing(lty) ) c(3, 3)             else rep(lty , length = 2)
        pch  <- if (missing(pch) ) c(21, 21)           else rep(pch , length = 2)
        col  <- if (missing(col) ) c("black", "red")   else rep(col , length = 2)
        bg   <- if (missing(bg)  ) c("white", "white") else rep(bg  , length = 2)

        x <- matrix(rep(object$size, each = object$nbest), nrow = object$nbest)
        y <- object$summary$aic[, object$size, drop = FALSE]

        par(mar = c(5, 4, 4, 4) + 0.1)

        plot.lmSubsets(object, main = main, sub = sub, xlab = xlab,
                       ylab = ylab[1], type = type[1], lty = lty[1],
                       pch = pch[1], col = col[1], bg = bg[1], ...,
                       legend = NULL)

        if (!is.null(legend)) {
            legend("topright", legend = legend, lty = lty,
                   pch = pch, col = col, pt.bg = bg)
        }

        xlim <- range(object$size)
        ylim <- range(y[is.finite(y)])
        plot.window(xlim = xlim, ylim = ylim)

        axis(side = 4, at = pretty(ylim))
        mtext(ylab[2], side = 4, line = 3)

        matplot(x = x, y = y, type = type[2], lty = lty[2], pch = pch[2],
                col = col[2], bg = bg[2], ..., add = TRUE)
    }

    if (missing(legend)) legend <- c("Deviance (RSS)", paste("Value (", x$summary$penalty, ")", sep = ""))

    localPlot(x, ...)

    invisible(x)
}


## plot 'lmSelect' summary
##
## Args:
##   x      - (mcsSubset)
##   ...    - ignored
##   legend - (character[])
##
## Rval: (summary.lmSelect) invisible
##
## All arguments are passed to 'plot.default'.
##
plot.summary.lmSelect <- function (x, ..., legend) {
    localPlot <- function (object, main, sub = NULL, xlab, ylab,
                           type, lty, pch, col, bg, ...) {
        n <- length(object$summary$penalty)

        if (missing(main)) main <- "All subsets (summary)"

        if (missing(type)) {
            type <- c("o", "o", rep("o", n))
        } else {
            type <- rep(type, length = n + 2)
        }
        if (missing(lty)) {
            lty <- c(3, 1, rep(1, n))
        } else {
            lty <- rep(lty, length = n + 2)
        }
        if (missing(pch)) {
            pch <- c(21, 21, rep(21, n))
        } else {
            pch <- rep(pch, length = n + 2)
        }
        if (missing(col)) {
            col <- c("black", "black", rep(c("red", "green3", "cyan", "blue", "magenta"), n))
        } else {
            col <- rep(col, length = n + 2)
        }
        if (missing(bg)) {
            bg <- c("white", "white", rep("white", n))
        } else {
            bg <- rep(bg, length = n + 2)
        }

        aic <- c(object$aic, object$summary$aic)
        ylim <- range(aic[is.finite(aic)])
        plot.lmSelect(object, main = main, sub = sub, xlab = xlab,
                      ylab = ylab, type = type[1:2], lty = lty[1:2], pch = pch[1:2],
                      col = col[1:2], bg = bg[1:2], ..., ylim2 = ylim, legend = NULL)

        if (!is.null(legend)) {
            legend("topleft", legend = legend, lty = lty,
                   pch = pch, col = col, pt.bg = bg)
        }

        y <- object$summary$aic
        matplot(t(y), type = type[-(1:2)], lty = lty[-(1:2)], pch = pch[-(1:2)],
                col = col[-(1:2)], bg = bg[-(1:2)], ..., add = TRUE)
    }

    if (missing(legend)) legend <- c("Deviance (RSS)", paste("Value (", x$penalty, ")", sep = ""),
                                     paste("Value (", x$summary$penalty, ", k = ", format(attr(x$summary$penalty, "k"), nsmall = 2), ")", sep = ""))

    localPlot(x, ...)

    invisible(x)
}
