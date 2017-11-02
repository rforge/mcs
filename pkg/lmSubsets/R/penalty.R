

penalty <- function (spec, ...) {
    if (is.function(spec)) {
        s <- evalq(substitute(penalty), parent.frame())
        s <- deparse(s)

        format <- function (...)  s
        eval <- function (obj, ...) {
            ans <- obj$rss
            ans[] <- mapply(spec, obj$rank, ans)
            ans
        }
        arg <- function (...)  spec
    } else if (tolower(spec) == "aic") {
        format <- function (...)  "AIC"
        eval <- function (obj, ...) {
            ll <- stats_log_lik(obj$nobs, obj$weights, obj$rss)
            stats_aic(ll, 2, obj$rank + 1)
        }
        arg <- function (...)  2
    } else if (tolower(spec) == "bic") {
        format <- function (...)  "BIC"
        eval <- function (obj, ...) {
            ll <- stats_log_lik(obj$nobs, obj$weights, obj$rss)
            stats_bic(ll, obj$nobs, obj$rank + 1)
        }
        arg <- function (nobs, ...)  log(nobs)
    } else if (is.numeric(spec)) {
        format <- function (...) {
            paste0("AIC (k = ", format_numeric(spec, ...), ")")
        }
        eval <- function (obj, ...) {
            ll <- stats_log_lik(obj$nobs, obj$weights, obj$rss)
            stats_aic(ll, spec, obj$rank + 1)
        }
        arg <- function (...)  spec
    } else {
        stop ("invalid penalty")
    }

    list(format = format, eval = eval, arg = arg)
}


format_penalty <- function (x, ...) {
    x$format(...)
}


eval_penalty <- function (x, ...) {
    x$eval(...)
}


arg_penalty <- function (x, ...) {
    x$arg(...)
}
