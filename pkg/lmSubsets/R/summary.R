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
    ## aic
    N <- object$nobs
    if (is.null(w <- object$weights)) {
        S <- 0
    } else {
        S <- sum(log(w))
    }

    pen <- NULL
    aic <- NULL
    for (p in penalty) {
        if (tolower(p) == "aic") {
            k <- 2
            pen <- structure(c(pen, "AIC"),
                             k = c(attr(pen, "k"), k))
        } else if (tolower(p) == "bic") {
            k <- log(N)
            pen <- structure(c(pen, "BIC"),
                             k = c(attr(pen, "k"), k))
        } else {
            k <- as.numeric(p)
            pen <- structure(c(pen, "AIC"),
                             k = c(attr(pen, "k"), k))
        }
    
        ll <- 0.5 * (S - N * (log(2 * pi) + 1 - log(N) + log(object$rss)))
        aic <- rbind(aic, -2 * ll + k * object$df)
    }

    object$summary <- list(penalty = pen,
                           aic = matrix(aic, ncol = object$nbest))

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
    fit <- format(x$summary$aic[, x$size, drop = FALSE], nsmall = 2)
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
##   ...    - ignored
##
## Rval: (summary.lmSubsets) invisible
##
## All arguments are passed to 'plot.default'.
##
plot.summary.lmSubsets <- function (x, ...) {
    ## start
    plot.new()
    box()

    ## title
    title(main = "lmSubsets (summary)")

    ## legend
    legend("topright", legend = c("Value", "Deviance"),
           lty = 3, pch = 21,
           col = c("red", "black"), pt.bg = "white")
    
    ## plot deviance
    xlim <- range(x$size)
    ylim <- range(x$rss[!is.na(x$rss)])
    plot.window(xlim = xlim, ylim = ylim)
    axis(1, at = pretty(xlim))
    axis(4, at = pretty(ylim))
    mtext("Size", side = 1, line = 3, col = "black")
    mtext("Deviance", side = 4, line = 3, col = "black")
    plot.lmSubsets(x, add = TRUE)

    ## plot value
    ylim <- range(x$summary$aic[!is.na(x$summary$aic)])
    plot.window(xlim = xlim, ylim = ylim)
    axis(2, at = pretty(ylim))
    mtext("Value", side = 2, line = 3, col = "black")
    for (i in x$size) {
        lines(x = rep(i, x$nbest), y = x$summary$aic[, i],
              type = "o", lty = 3, pch = 21,
              col = "red", bg = "white")
    }


    ## done
    invisible(x)
}


## plot 'lmSelect' summary
##
## Args:
##   x      - (mcsSubset)
##   ...    - ignored
##
## Rval: (summary.lmSelect) invisible
##
## All arguments are passed to 'plot.default'.
##
plot.summary.lmSelect <- function (x, ...) {
    par(mar = c(5, 4, 4, 4) + 0.1)

    ## start
    plot.new()
    box()

    xlim <- c(1, x$nbest)

    ## title
    title(main = "lmSelect (summary)")

    ## legend
    legend("topleft", legend = c("Value", "Deviance"),
           lty = c(1, 3), pch = 21,
           col = c("red", "black"), pt.bg = c("red", "white"),
           text.col = c("red", "black"))

    ## plot deviance
    plot.window(xlim = xlim, ylim = range(x$rss))
    lines(x$rss, type = "o", lty = 3, pch = 21, col = "black", bg = "white")
    axis(4, at = pretty(range(x$rss)))
    mtext("Deviance", side = 4, line = 3, col = "black")

    ## plot value
    plot.window(xlim = xlim, ylim = range(x$summary$aic))
    matplot(t(x$summary$aic), type = "o", lty = 1, pch = 21, col = 2:7, bg = 2:7, add = TRUE)
    axis(2, at = pretty(range(x$summary$aic)))
    mtext("Value", side = 2, line = 3)

    ## x axis
    axis(1, at = pretty(xlim))
    mtext("Best", side = 1, line = 3)

    ## done
    invisible(x)
}
