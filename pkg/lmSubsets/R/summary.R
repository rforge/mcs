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
print.summary.lmSubsets <- function (x, ...) {
    catln <- function (..., sep = "") base::cat(..., "\n", sep = sep)
    paste <- function (..., sep = "") base::paste(..., sep = sep)

    ## call
    catln()
    catln("Call:")
    catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))

    ## summary
    catln()
    cat("Summary:")
    output <- format(attr(x$summary$penalty, "k"), nsmall = 2)
    output <- paste(x$summary$penalty, " (penalty = ", output, ")")
    output <- as.matrix(output, ncol = 1)
    colnames(output) <- ""
    rownames(output) <- "  Value:"
    print(output, quote = FALSE)

    ## fit
    catln()
    catln("Model fit:")
    catln("  best x size")
    catln("    DEVIANCE")
    catln("    Summary")
    output <- x$rss[, x$nmin:x$nmax, drop = FALSE]
    output <- rbind(output, x$summary$val[, x$nmin:x$nmax, drop = FALSE])
    output <- output[order(c(1:x$nbest, 1:x$nbest)), ]
    row.dev <- seq(1, by = 2, length.out = x$nbest)
    row.sum <- seq(2, by = 2, length.out = x$nbest)
    star <- apply(output[row.sum, , drop = FALSE], 1, which.min)
    star <- cbind(row.sum, star)
    output <- ifelse(is.na(output), "", format(output, nsmall = 2))
    output[star] <- paste(output[star], "*", sep = "")
    rownames(output)[row.dev] <- paste("  ", rownames(output)[row.dev])
    rownames(output)[row.sum] <- rep("", x$nbest)
    print(output, quote = FALSE)

    ## pad
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

    N <- length(x$summary$penalty)

    ## call
    catln()
    catln("Call:")
    catln("  ", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))

    ## summary
    catln()
    cat("Summary:")
    output <- format(attr(x$summary$penalty, "k"), nsmall = 2)
    output <- paste("(penalty = ", output, ")")
    output <- paste(x$summary$penalty, output, sep = " ")
    output <- as.matrix(output, ncol = 1)
    colnames(output) <- ""
    rownames(output) <- c("  Value:", rep("", N - 1))
    print(output, quote = FALSE)

    ## fit
    catln()
    catln("Model fit:")
    output <- rbind(x$df, x$rss, x$val, x$summary$val)
    star <- apply(output[-1, , drop = FALSE], 1, which.min)
    star <- cbind(2:nrow(output), star)
    output <- format(output, nsmall = 2)
    output[star] <- paste(output[star], "*", sep = "")
    rownames(output) <- paste("  ", c("df", "Deviance", "VALUE", "Summary", rep("", N - 1)))
    print(output, quote = FALSE)

    ## pad
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
    localPlot <- function (object, main, sub, xlab, ylab,
                           type.sum = "o", lty.sum = c(1, 3),
                           pch.sum = c(16, 21), col.sum = "red",
                           bg.sum = "white", type = "o", lty = c(1, 3),
                           pch = c(16, 21), col = "black", bg = "white", ...) {
        if (missing(main))  main <- "All subsets (summary)"
        if (missing(xlab))  xlab <- "Number of regressors"
        if (missing(ylab))  ylab <- c("Deviance", "Value")

        type.sum <- rep(type.sum, length = 2)
        lty.sum  <- rep(lty.sum , length = 2)
        pch.sum  <- rep(pch.sum , length = 2)
        col.sum  <- rep(col.sum , length = 2)
        bg.sum   <- rep(bg.sum  , length = 2)

        x <- matrix(rep(object$nmin:object$nmax, each = object$nbest),
                    nrow = object$nbest)
        y <- object$summary$val[object$nbest:1, object$nmin:object$nmax,
                                drop = FALSE]

        par(mar = c(5, 4, 4, 4) + 0.1)

        plot.lmSubsets(object, main = main, sub = sub, xlab = xlab, ylab = ylab[1],
                       type = type, lty = lty, pch = pch,
                       col = col, bg = bg, ..., legend = NULL)

        if (!is.null(legend)) {
            legend("topright", legend = legend, lty = c(lty[1], lty.sum[1]),
                   pch = c(pch[1], pch.sum[1]), col = c(col[1], col.sum[1]),
                   pt.bg = c(bg[1], bg.sum[1]), bty = "n")
        }

        xlim <- c(object$nmin, object$nmax)
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

    ## default legend
    if (missing(legend)) {
        legend <- attr(x$summary$penalty, "k")
        legend <- format(legend, nsmall = 2)
        legend <- paste("Summary (", x$summary$penalty, ", penalty = ",
                        legend, ")", sep = "")
        legend <- c("Deviance (RSS)", legend)
    }

    ## plot
    localPlot(x, ...)

    ## done
    invisible(x)
}


## plot 'lmSelect' summary
##
## Args:
##   x      - (mcsSubset)
##   ...    - forwarded to plotting functions
##   legend - (character[])
##
## Rval: (summary.lmSelect) invisible
##
## All arguments are passed to 'plot.default'.
##
plot.summary.lmSelect <- function (x, ..., legend) {
    localPlot <- function (object, main, sub, xlab, ylab, type.sum = "o",
                           lty.sum = 1, pch.sum = 21, col.sum, bg.sum = "white",
                           type = c("o", "o"), lty = c(3, 1), pch = c(21, 16),
                           col = c("black", "black"), bg = c("white", "white"),
                           ...) {
        n <- length(object$summary$penalty)

        if (missing(main))  main <- "All subsets (summary)"
        if (missing(xlab))  xlab <- "Number of regressors"
        if (missing(ylab))  ylab <- c("Deviance", "Value")

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

    ## default legend
    if (missing(legend)) {
        legend.2 <- format(attr(x$penalty, "k"), nsmall = 2)
        legend.2 <- paste("Value (", x$penalty, ", penalty = ",
                          legend.2, ")", sep = "")
        legend.3 <- format(attr(x$summary$penalty, "k"), nsmall = 2)
        legend.3 <- paste("Summary (", x$summary$penalty, ", penalty = ",
                          legend.3, ")", sep = "")
        legend <- c("Deviance (RSS)", legend.2, legend.3)
    }

    ## plot
    localPlot(x, ...)

    ## done
    invisible(x)
}
