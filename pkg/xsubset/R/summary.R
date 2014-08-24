

#############
## SUMMARY ##
#############


## summary for 'xsubset' objects
##
## Args:
##   object  - (xsubset)
##   rank    - (integer)
##   penalty - (numeric)
##   ...     - ignored
##
## Rval: (summary.xsubset)
##
summary.xsubset <- function (object, penalty = 2, ...) {
  paste <- function (..., sep = "") base::paste(..., sep = sep)

  ## penalty
  object$penalty <- penalty
  ## aic
  aic <- AIC(object, size = 1:object$nvar, rank = 1:object$nbest,
             k = penalty)
  aic <- matrix(aic$AIC, nrow = object$nbest, ncol = object$nvar)

  ## order
  pi <- order(aic)
  pi <- array(c((pi - 1) %%  object$nbest + 1,
                (pi - 1) %/% object$nbest + 1),
              dim = c(length(aic), 2))
  nok <- is.na(aic[pi])
  pi <- pi[!nok, , drop = FALSE]
  ## nbest
  object$nbest <- nrow(pi)

  cnames <- paste(1:object$nbest, ".")
  ## tolerance
  object$tolerance <- array(object$tolerance[pi[, 2]], dim = object$nbest,
                            dimnames = list(cnames))
  ## rss
  object$rss <- array(object$rss[pi], dim = object$nbest,
                      dimnames = list(cnames))
  ## aic
  object$aic <- array(aic[pi], dim = object$nbest,
                      dimnames = list(cnames))
  ## which table
  sel <- cbind(rep(1:object$nvar, nrow(pi)),
               rep(pi[, 1], each = object$nvar),
               rep(pi[, 2], each = object$nvar))
  object$which <- array(object$which[sel], dim = c(object$nvar, object$nbest),
                        dimnames = list(object$x.names, cnames))
  ## size
  object$size <- array(apply(object$which, 2, sum), dim = object$nbest,
                       dimnames = list(cnames))

  ## class
  class(object) <- c("summary.xsubset", "xsubset2", "xsubset0")

  ## done
  object
}


## summary for 'xsubset2' objects
##
## Args:
##   object  - (xsubset2)
##   penalty - (numeric)
##   ...     - ignored
##
## Rval:  (summary.xsubset2)
##
## Do nothing for objects of type 'xsubset2'.
##
summary.xsubset2 <- function (object, penalty = NULL,...) {
  paste <- function (..., sep = "") base::paste(..., sep = sep)

  if (is.null(penalty)) {
    class(object) <- c("summary.xsubset", "xsubset2", "xsubset0")
    return (object)
  }

  ## penalty
  object$penalty <- penalty
  ## aic
  aic <- AIC(object, rank = 1:object$nbest, k = penalty)
  aic <- aic$AIC

  ## order
  pi <- order(aic)
  nok <- is.na(aic[pi])
  pi <- pi[!nok]
  ## nbest
  object$nbest <- length(pi)

  cnames <- paste(1:object$nbest, ".")
  ## tolerance
  object$tolerance <- array(object$tolerance[pi], dim = object$nbest,
                            dimnames = list(cnames))
  ## rss
  object$rss <- array(object$rss[pi], dim = object$nbest,
                      dimnames = list(cnames))
  ## aic
  object$aic <- array(aic[pi], dim = object$nbest,
                      dimnames = list(cnames))
  ## which table
  object$which <- array(object$which[, pi], dim = c(object$nvar, object$nbest),
                        dimnames = list(object$x.names, cnames))
  ## size
  object$size <- array(apply(object$which, 2, sum), dim = object$nbest,
                       dimnames = list(cnames))

  ## class
  class(object) <- c("summary.xsubset", "xsubset2", "xsubset0")

  ## done
  object
}


## print 'xsubset' summary
##
## Args:
##   x      - (summary.xsubset)
##   digits - (integer)
##   ...    - ignored
##
## Rval: (summary.xsubset) invisible
##
print.summary.xsubset <- function (x, digits = NULL, ...)
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
  catln("Selected variables (by rank):")
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
    
  invisible(x)
}


## plot 'xsubset' summary
##
## Args:
##   x      - (xsubset)
##   type   - (character) plot type
##   main   - (character) main title
##   xlab   - (character) x label
##   ylab   - (character) y label
##   col    - (integer[]|character[]) color
##   lty    - (integer) line type
##   legend - (logical) legend?
##   ...    - forwarded
##
## Rval: (summary.xsubset) invisible
##
## All arguments are passed to 'plot.default'.
##
plot.summary.xsubset <- function (x, type = "b", main = NULL, xlab = NULL,
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
    xlab <- "Rank"
  }

  ## color
  col <- rep(col, length.out = 2)

  ## line type
  lty <- rep(lty, length.out = 2)

  ## rank
  rank <- 1:x$nbest

  ## plot rss
  plot(rank, x$rss, type = type[1], main = main, sub = sub,
       xlab = xlab, ylab = ylab, col = col[1], lty = lty[1], ...)

  ## plot aic
  new.old <- getOption("new")
  par(new = TRUE)
  plot(rank, x$aic, type = type[2], xlab = "", ylab = "",
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
