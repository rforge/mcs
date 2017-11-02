##
## File:  summary.R
##



#################
##  LMSUBSETS  ##
#################


## summary for 'lmSubsets' objects
##
3## Arguments:
##   object - (lmSubsets)
##   ...    - ignored
##
## Result: (summary.lmSubsets)
##
summary.lmSubsets <- function (object, ...) {
    s <- stats_lmSubsets(object)

    ans <- object[c("call", "terms", "nvar", "nbest", "sizes",
                    if (!is.null(object$weights)) "weights")]

    ans$sigma <- sigma_stats(s)
    ans$r2 <- r2_stats(s)
    ans$r2_adj <- r2_adj_stats(s, r2 = ans$r2)
    ans$p_value <- pval_stats(s)
    ans$mallows_cp <- cp_stats(s)
    ans$aic <- aic_stats(s)
    ans$bic <- bic_stats(s)

    class(ans) <- "summary.lmSubsets"

    ans
}


## format 'lmSubsets' summary
##
## Arguments:
##   x   - (summary.lmSubsets)
##   ... - forwarded
##
## Result: (list)
##
format.summary.lmSubsets <- function (x, ...) {
    ans <- list()

    ans$call <- paste(deparse(x$call), sep = "\n", collapse = "\n")

    ans$sigma <- x$sigma[, x$sizes, drop = FALSE]
    ans$sigma <- ifelse(is.na(ans$sigma), "", format_numeric(ans$sigma, ...))

    ans$r2 <- rbind(x$r2, x$r2_adj)
    ans$r2 <- ans$r2[order(rep(seq_len(x$nbest), 2)), x$sizes, drop = FALSE]
    ans$r2 <- ifelse(is.na(ans$r2), "", format_numeric(ans$r2, ...))
    rownames(ans$r2)[2 * seq_len(x$nbest)] <- ""

    ans$p_value <- x$p_value[, x$sizes, drop = FALSE]
    ans$p_value <- ifelse(is.na(ans$p_value), "", format.pval(ans$p_value))

    ans$mallows_cp <- x$mallows_cp[, x$sizes, drop = FALSE]
    ans$mallows_cp <- ifelse(is.na(ans$mallows_cp), "", format_numeric(ans$mallows_cp, ...))

    ans$aic <- rbind(x$aic, x$bic)
    ans$aic <- ans$aic[order(rep(seq_len(x$nbest), 2)), x$sizes, drop = FALSE]
    ans$aic <- ifelse(is.na(ans$aic), "", format_numeric(ans$aic, ...))
    rownames(ans$aic)[2 * seq_len(x$nbest)] <- ""

    ans
}


## print 'lmSubsets' summary
##
## Arguments:
##   x   - (summary.lmSubsets)
##   ... - forwarded
##
## Result: (summary.lmSubsets) invisible
##
print.summary.lmSubsets <- function (x, ...) {
    catln <- function (...)  base::cat(..., "\n", sep = "")
    print <- function (..., skip = 0, indent = 0) {
        output <- capture.output(base::print(...))
        if (skip > 0)  output <- output[-seq_len(skip)]
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

    ## sigma
    catln("Residual standard deviation:")
    catln("  best x size")
    catln("    sigma")
    print(fmt$sigma, quote = FALSE, indent = 2)

    catln()

    ## R squared
    catln("Coefficient of determination:")
    catln("  best x size")
    catln("    R^2")
    catln("    R^2 adj.")
    print(fmt$r2, quote = FALSE, indent = 2)

    catln()

    ## p-value
    catln("P-value (F-statistic):")
    catln("  best x size")
    catln("    p-value")
    print(fmt$p_value, quote = FALSE, indent = 2)

    catln()
    
    ## mallows' cp
    catln("Mallows' Cp:")
    catln("  best x size")
    catln("    Cp")
    print(fmt$mallows_cp, quote = FALSE, indent = 2)

    catln()

    ## aic
    catln("Akaike's information criterion:")
    catln("  best x size")
    catln("    AIC")
    catln("    BIC")
    print(fmt$aic, quote = FALSE, indent = 2)

    catln()

    ## done
    invisible(x)
}



################
##  LMSELECT  ##
################


## summary for 'lmSelect' objects
##
## Arguments:
##   object  - (lmSelect)
##   ...     - ignored
##
## Result: (summary.lmSelect)
##
summary.lmSelect <- function (object, ...) {
    s <- stats_lmSelect(object)

    ans <- object[c("call", "terms", "nvar", "nbest", "sizes",
                    if (!is.null(object$weights)) "weights")]

    ans$rank <- s$rank
    ans$sigma <- sigma_stats(s)
    ans$r2 <- r2_stats(s)
    ans$r2_adj <- r2_adj_stats(s, r2 = ans$r2)
    ans$p_value <- pval_stats(s)
    ans$mallows_cp <- cp_stats(s)
    ans$aic <- aic_stats(s)
    ans$bic <- bic_stats(s)

    class(ans) <- "summary.lmSelect"

    ans
}


## format 'lmSelect' summary
##
## Arguments:
##   x   - (summary.lmSelect)
##   ... - forwarded
##
## Result: (list)
##
format.summary.lmSelect <- function (x, ...) {
    ans <- list()

    ans$call <- paste(deparse(x$call), sep = "\n", collapse = "\n")

    names <- format_ordinal(seq_len(x$nbest))

    ans$rank <- t(x$rank)
    rownames(ans$rank) <- "rank"
    colnames(ans$rank) <- names

    ans$sigma <- t(x$sigma)
    ans$sigma <- format_numeric(ans$sigma, ...)
    rownames(ans$sigma) <- "sigma"
    colnames(ans$sigma) <- names

    ans$r2 <- rbind(x$r2, x$r2_adj)
    ans$r2 <- format_numeric(ans$r2, ...)
    rownames(ans$r2) <- c("R2", "R2 adj.")
    colnames(ans$r2) <- names

    ans$p_value <- format.pval(x$p_value)
    ans$p_value <- t(ans$p_value)
    rownames(ans$p_value) <- "p-val."
    colnames(ans$p_value) <- names

    ans$mallows_cp <- t(x$mallows_cp)
    ans$mallows_cp <- format_numeric(ans$mallows_cp, ...)
    rownames(ans$mallows_cp) <- "Cp"
    colnames(ans$mallows_cp) <- names

    ans$aic <- rbind(x$aic, x$bic)
    ans$aic <- format_numeric(ans$aic, ...)
    rownames(ans$aic) <- c("AIC", "BIC")
    colnames(ans$aic) <- names

    ans
}


## print 'lmSelect' summary
##
## Arguments:
##   x   - (summary.lmSelect)
##   ... - forwarded
##
## Result: (summary.lmSelect) invisible
##
print.summary.lmSelect <- function (x, ...) {
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

    ## rank
    catln("Rank:")
    catln("  best")
    print(fmt$rank, quote = FALSE, indent = 2)

    catln()

    ## sigma
    catln("Residual standard deviation:")
    catln("  best")
    print(fmt$sigma, quote = FALSE, indent = 2)

    catln()

    ## R squared
    catln("Coefficient of determination:")
    catln("  best")
    print(fmt$r2, quote = FALSE, indent = 2)

    catln()
    
    ## mallows' cp
    catln("Mallows' Cp:")
    catln("  best")
    print(fmt$mallows_cp, quote = FALSE, indent = 2)

    catln()

    ## p-value
    catln("P-value (F-statistic):")
    catln("  best")
    print(fmt$p_value, quote = FALSE, indent = 2)

    catln()

    ## aic
    catln("Akaike's information criterion:")
    catln("  best")
    print(fmt$aic, quote = FALSE, indent = 2)

    catln()

    ## done
    invisible(x)
}
