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
summary.lmSubsets <- function (object, penalty = "BIC", ...) {

    ## aic
    N <- object$nobs
    if (is.null(w <- object$weights)) {
        S <- 0
    } else {
        S <- sum(log(w[w != 0]))
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

    object$summary <- list(penalty = penalty, val = aic)

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
summary.lmSelect <- function (object, penalty = "AIC", ...) {
    ## log lik
    N <- object$nobs
    if (is.null(w <- object$weights)) {
        S <- 0
    } else {
        S <- sum(log(w[w != 0]))
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
                           val     = aic                  )

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
    val <- x$summary$val[, x$size, drop = FALSE]
    fit <- ifelse(is.na(val), "", format(val, nsmall = 2))
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
    fit <- format(rbind(x$df, x$rss, x$summary$val), nsmall = 2)
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
##   ...    - passed to plotting functions
##   legend - (character[])
##
## Rval: (summary.lmSubsets) invisible
##
plot.summary.lmSubsets <- function (x, ..., legend) {
    localPlot <- function (object, main, sub = NULL, xlab, ylab,
                           type.sum = "o", lty.sum = c(1, 3),
                           pch.sum = c(16, 21), col.sum = "red",
                           bg.sum = "white", type = "o", lty = c(1, 3),
                           pch = c(16, 21), col = "black", bg = "white", ...) {
        if (missing(main)) main <- "All subsets (summary)"
        if (missing(xlab)) xlab <- "Number of regressors"
        if (missing(ylab)) ylab <- c("Deviance", "Value")

        type.sum <- rep(type.sum, length = 2)
        lty.sum  <- rep(lty.sum , length = 2)
        pch.sum  <- rep(pch.sum , length = 2)
        col.sum  <- rep(col.sum , length = 2)
        bg.sum   <- rep(bg.sum  , length = 2)

        x <- matrix(rep(object$size, each = object$nbest), nrow = object$nbest)
        y <- object$summary$val[object$nbest:1, object$size, drop = FALSE]

        par(mar = c(5, 4, 4, 4) + 0.1)

        plot.lmSubsets(object, main = main, sub = sub, xlab = xlab, ylab = ylab,
                       type = type, lty = lty, pch = pch,
                       col = col, bg = bg, ..., legend = NULL)

        if (!is.null(legend)) {
            legend("topright", legend = legend, lty = c(lty[1], lty.sum[1]),
                   pch = c(pch[1], pch.sum[1]), col = c(col[1], col.sum[1]),
                   pt.bg = c(bg[1], bg.sum[1]), bty = "n")
        }

        xlim <- range(object$size)
        ylim <- range(y[is.finite(y)])
        plot.window(xlim = xlim, ylim = ylim)

        axis(side = 4, at = pretty(ylim))
        mtext(ylab[2], side = 4, line = 3)

        matplot(x = x, y = y, type = type.sum[2], lty = lty.sum[2],
                pch = pch.sum[2], col = col.sum[2], bg = bg.sum[2], ...,
                add = TRUE)
        lines(x[1, ], y[object$nbest, ], type = type.sum[1], lty = lty.sum[1],
              pch = pch.sum[1], col = col.sum[1], bg = bg.sum[1], ...)
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
    localPlot <- function (object, main, sub = NULL, xlab, ylab, type.sum = "o",
                           lty.sum = 1, pch.sum = 21, col.sum, bg.sum = "white",
                           type = c("o", "o"), lty = c(3, 1), pch = c(21, 16),
                           col = c("black", "black"), bg = c("white", "white"),
                           ...) {
        n <- length(object$summary$penalty)

        if (missing(main)) main <- "All subsets (summary)"
        if (missing(xlab)) xlab <- "Number of regressors"
        if (missing(ylab)) ylab <- c("Deviance", "Value")

        if (missing(col.sum)) {
            col.sum <- c("red", "green3", "cyan", "blue", "magenta")
        }

        type.sum <- rep(type.sum, length = n)
        lty.sum  <- rep(lty.sum , length = n)
        pch.sum  <- rep(pch.sum , length = n)
        col.sum  <- rep(col.sum , length = n)
        bg.sum   <- rep(bg.sum  , length = n)

        val <- c(object$val, object$summary$val)
        ylim <- range(val[is.finite(val)])
        plot.lmSelect(object, main = main, sub = sub, xlab = xlab, ylab = ylab,
                      type = type, lty = lty, pch = pch, col = col, bg = bg,
                      ..., ylim2 = ylim, legend = NULL)

        if (!is.null(legend)) {
            legend("topleft", legend = legend, lty = c(lty, lty.sum),
                   pch = c(pch, pch.sum), col = c(col, col.sum),
                   pt.bg = c(bg, bg.sum), bty = "n")
        }

        y <- object$summary$val
        matplot(t(y), type = type.sum, lty = lty.sum, pch = pch.sum,
                col = col.sum, bg = bg.sum, ..., add = TRUE)
    }

    if (missing(legend)) legend <- c("Deviance (RSS)",
                                     paste("Value (", x$penalty, ", k = ", format(attr(x$penalty, "k"), nsmall = 2), ")", sep = ""),
                                     paste("Value (", x$summary$penalty, ", k = ", format(attr(x$summary$penalty, "k"), nsmall = 2), ")", sep = ""))

    localPlot(x, ...)

    invisible(x)
}
