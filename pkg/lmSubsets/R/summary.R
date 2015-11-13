

#############
## SUMMARY ##
#############


## summary for 'lmSubsets' objects
##
## Args:
##   object  - (lmSubsets)
##   penalty - (numeric|"AIC"|"BIC") passed to 'AIC.lmSubsets'
##   ...     - ignored
##
## Rval: (summary.lmSubsets)
##
summary.lmSubsets <- function (object, penalty = 2, ...) {
  paste <- function (..., sep = "") base::paste(..., sep = sep)

  x.names <- variable.names(object, .full = TRUE)

  ## aic
  aic <- AIC(object, size = 1:object$nvar, best = 1:object$nbest,
             k = penalty)
  aic <- matrix(aic$AIC, nrow = object$nbest)
  ## penalty
  object$penalty <- penalty
  ## order
  pi <- order(aic)
  pi <- array(c((pi - 1) %%  object$nbest + 1,
                (pi - 1) %/% object$nbest + 1),
              dim = c(length(aic), 2))
  nok <- is.na(aic[pi])
  pi <- pi[!nok, , drop = FALSE]
  ## nbest
  object$nbest <- nrow(pi)

  .cnames <- paste(1:object$nbest, ".")
  ## tolerance
  object$tolerance <- array(object$tolerance[pi[, 2]], dim = object$nbest,
                            dimnames = list(.cnames))
  ## rss
  object$rss <- array(object$rss[pi], dim = object$nbest,
                      dimnames = list(.cnames))
  ## aic
  object$aic <- array(aic[pi], dim = object$nbest,
                      dimnames = list(.cnames))
  ## which table
  sel <- cbind(rep(1:object$nvar, nrow(pi)),
               rep(pi[, 1], each = object$nvar),
               rep(pi[, 2], each = object$nvar))
  object$which <- array(object$which[sel], dim = c(object$nvar, object$nbest),
                        dimnames = list(x.names, .cnames))
  ## size
  object$size <- array(apply(object$which, 2, sum), dim = object$nbest,
                       dimnames = list(.cnames))

  ## class
  class(object) <- "summary.lmSubsets"

  ## done
  object
}


## print 'lmSubsets' summary
##
## Args:
##   x      - (summary.lmSubsets)
##   digits - (integer)
##   ...    - ignored
##
## Rval: (summary.lmSubsets) invisible
##
print.summary.lmSubsets <- function (x, digits = NULL, ...)
{
    catln <- function (..., sep = "") base::cat(..., "\n", sep = sep)
    paste <- function (..., sep = "") base::paste(..., sep = sep)

    ## digits
    if (is.null(digits)) {
        digits <- max(3, getOption("digits") - 3)
    }
    ## call
    catln()
    catln("Call:")
    catln(deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)))
    ## variable table
    catln()
    catln("Selected variables (best first):")
    which.x <- ifelse(x$which, "x", "")
    rownames(which.x)[x$include] <- paste("+", rownames(which.x)[x$include])
    rownames(which.x)[x$exclude] <- paste("-", rownames(which.x)[x$exclude])
    print(which.x, quote = FALSE)
    ## fit
    catln()
    catln("Model fit:")
    fit <- format(rbind(x$aic, x$rss), digits = digits)
    fit <- rbind(fit, format(x$size))
    rownames(fit) <- c("AIC", "RSS", "(size)")
    print(fit, quote = FALSE)
    catln()
    catln("AIC: k = ", format(x$penalty, digits = digits))

    ## done
    invisible(x)
}


## plot 'lmSubsets' summary
##
## Args:
##   x      - (lmSubsets)
##   type   - (character) plot type
##   main   - (character) main title
##   xlab   - (character) x label
##   ylab   - (character) y label
##   col    - (integer[]|character[]) color
##   lty    - (integer) line type
##   legend - (logical) legend?
##   ...    - forwarded
##
## Rval: (summary.lmSubsets) invisible
##
## All arguments are passed to 'plot.default'.
##
plot.summary.lmSubsets <- function (x, type = "b", main = NULL, xlab = NULL,
                                    ylab = "", col = c("blue", "red"), lty = 1,
                                    legend = TRUE, ...)
{
    digits <- max(3, getOption("digits") - 3)

    ## type
    type <- rep(type, length.out = 2)
    ## main title
    if (is.null(main)) {
        main <- paste("AIC and residual sum of squares")
    }
    ## sub title
    sub <- paste("AIC (k = ", format(x$penalty, digits = digits),
                 ")", sep = "")
    ## x label
    if (is.null(xlab)) {
        xlab <- "Best"
    }
    ## color
    col <- rep(col, length.out = 2)
    ## line type
    lty <- rep(lty, length.out = 2)
    ## best
    best <- 1:x$nbest
    ## plot rss
    plot(best, x$rss, type = type[1], main = main, sub = sub,
         xlab = xlab, ylab = ylab, col = col[1], lty = lty[1], ...)
    ## plot aic
    new.old <- getOption("new")
    par(new = TRUE)
    plot(best, x$aic, type = type[2], xlab = "", ylab = "",
         col = col[2], lty = lty[2], axes = FALSE, ...)
    par(new = new.old)
    ## legend
    if (legend) {
        legend("top", c("RSS", "AIC"), lty = lty, col = col, bty = "n")
    }
    ## axes
    axis(4)

    ## done
    invisible(x)
}
