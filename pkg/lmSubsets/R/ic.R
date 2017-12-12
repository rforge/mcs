

ic <- function (penalty, ...) {
    if (is.function(penalty)) {
        f <- penalty

        s <- evalq(substitute(penalty), parent.frame())
        s <- deparse(s)

        format <- function (...)  s
        eval <- function (obj, ...) {
            rss <- deviance(obj, size = NULL, best = NULL, na.rm = FALSE)
            size <- obj$submodel$SIZE

            with(obj$submodel, mapply(f, size, rss))
        }
        penalty <- function (nobs, ...)  f
    } else if (tolower(penalty) == "aic") {
        format <- function (...)  "AIC"
        eval <- function (obj, ...) {
            rss <- deviance(obj, size = NULL, best = NULL, na.rm = FALSE)
            size <- obj$submodel$SIZE

            ll <- stats_log_lik(obj$nobs, obj$weights, rss)
            stats_aic(ll, 2, size + 1)
        }
        penalty <- function (nobs, ...)  2.0
    } else if (tolower(penalty) == "bic") {
        format <- function (...)  "BIC"
        eval <- function (obj, ...) {
            rss <- deviance(obj, size = NULL, best = NULL, na.rm = FALSE)
            size <- obj$submodel$SIZE

            ll <- stats_log_lik(obj$nobs, obj$weights, rss)
            stats_bic(ll, obj$nobs, size + 1)
        }
        penalty <- function (nobs, ...)  log(nobs)
    } else if (is.numeric(penalty)) {
        k <- penalty

        format <- function (...) {
            paste0("AIC (k = ", format_default(k, ...), ")")
        }
        eval <- function (obj, ...) {
            rss <- deviance(obj, size = NULL, best = NULL, na.rm = FALSE)
            size <- obj$submodel$SIZE

            ll <- stats_log_lik(obj$nobs, obj$weights, rss)
            stats_aic(ll, k, size + 1)
        }
        penalty <- function (nobs, ...)  k
    } else {
        stop ("invalid penalty")
    }

    list(format = format, eval = eval, penalty = penalty)
}


format_ic <- function (x, ...) {
    x$format(...)
}


eval_ic <- function (x, ...) {
    x$eval(...)
}


penalty_ic <- function (x, nobs, ...) {
    x$penalty(nobs, ...)
}
