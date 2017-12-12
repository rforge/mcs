

stats_formula <- function (x_names, y_name, intercept, env = NULL) {
    f <- paste0(x_names, collapse = " + ")
    f <- paste0(y_name, " ~ ", f)

    if (intercept) {
        f <- paste0(f, " + 1")
    } else {
        f <- paste0(f, " + 0")
    }

    if (is.null(env)) {
        env <- parent.frame(2L)
    }

    formula(f, env)
}


stats_coef <- function (x, y, o, w) {
    if (!is.null(o)) {
        y <- y - o
    }

    if (!is.null(w)) {
        ok <- w != 0

        if (any(!ok)) {
            w <- w[ok]
            x <- x[ok, , drop = FALSE]
            y <- y[ok]
        }

        w <- sqrt(w)
        x <- x * w
        y <- y * w
    }

    qr.solve(x, y)
}


stats_vcov <- function (x, y, o, w) {
    if (!is.null(o)) {
        y <- y - o
    }

    if (!is.null(w)) {
        ok <- w != 0

        if (any(!ok)) {
            w <- w[ok]
            x <- x[ok, , drop = FALSE]
            y <- y[ok]
        }

        w <- sqrt(w)
        x <- x * w
        y <- y * w
    }

    m <- nrow(x)
    n <- ncol(x)
    nn <- seq_len(n)

    xy <- cbind(x, y)
    qrz <- qr(xy)
    r <- qrz$qr[nn, nn, drop = FALSE]
    rdf <- m - n
    rss <- qrz$qr[n + 1, n + 1]^2
    sigma2 <- rss / rdf
    chol2inv(r) * sigma2
}


stats_fitted <- function (x, y, o, w) {
    if (!is.null(o)) {
        y <- y - o
    }

    if (!is.null(w)) {
        ok <- w != 0

        if (any(!ok)) {
            x0 <- x[!ok, , drop = FALSE]
            w <- w[ok]
            x <- x[ok, , drop = FALSE]
            y <- y[ok]
        }

        w <- sqrt(w)
        x <- x * w
        y <- y * w
    }

    qr <- qr(x)
    f <- qr.fitted(qr, y)

    if (!is.null(w)) {
        f <- f / w

        if (any(!ok)) {
            b <- qr.coef(qr, y)
            f0 <- drop(x0 %*% b)
            f <- c(f, f0)[order(c(which(ok), which(!ok)))]
        }
    }

    if (!is.null(o)) {
        f <- f + o
    }

    f
}


stats_residuals <- function (x, y, o, w) {
    if (!is.null(o)) {
        y <- y - o
    }

    if (!is.null(w)) {
        ok <- w != 0

        if (any(!ok)) {
            x0 <- x[!ok, , drop = FALSE]
            y0 <- y[!ok]
            w <- w[ok]
            x <- x[ok, , drop = FALSE]
            y <- y[ok]
        }

        w <- sqrt(w)
        x <- x * w
        y <- y * w
    }

    qr <- qr(x)
    r <- qr.resid(qr, y)

    if (!is.null(w)) {
        r <- r / w

        if (any(!ok)) {
            b <- qr.coef(qr, y)
            f0 <- drop(x0 %*% b)
            r0 <- y0 - f0
            r <- c(r, r0)[order(c(which(ok), which(!ok)))]
        }
    }

    r
}


stats_log_lik <- function (m, w, rss) {
    if (is.null(w)) {
        w <- rep.int(1, m)
    } else {
        w <- w[w != 0]
    }

    0.5 * (sum(log(w)) - m * (log(2 * pi) + 1 - log(m) + log(rss)))
}


stats_aic <- function (ll, k, df) {
    -2 * ll + k * df
}


stats_bic <- function (ll, m, df) {
    stats_aic(ll, log(m), df)
}


stats_mse <- function (m, rss) {
    rss / m
}


stats_sigma2 <- function (rss, rdf) {
    rss / rdf
}


stats_mss <- function (m, w, f, idf) {
    if (is.null(w)) {
        w <- rep.int(1, m)
    }

    if (idf > 0) {
        mu <- sum(w * f / sum(w))
    } else {
        mu <- 0
    }

    sum(w * (f - mu)^2)
}


stats_r2 <- function (rss, mss) {
    mss / (mss + rss)
}


stats_r2adj <- function (m, r2, idf, rdf) {
    1 - (1 - r2) * ((m - idf) / rdf)
}


stats_pval <- function (p, rss, mss, idf, rdf) {
    pf((mss / (p - idf)) / (rss / rdf), p - idf, rdf,
       lower.tail = FALSE)
}


stats_cp <- function (m, p, rss, sigma2_full) {
    rss / sigma2_full - m + 2 * p
}


stats_lmSubsets <- function (object, ...) {
    m <- object$nobs
    n <- object$nvar

    lm_full <- lm(object)
    rss_full <- deviance(lm_full)
    sigma2_full <- stats_sigma2(rss_full, m - n)

    idf <- if (object$intercept) 1 else 0

    o <- object$offset
    w <- object$weights

    p <- object$submodel$SIZE
    best <- object$submodel$BEST
    rss <- object$submodel$RSS

    y <- model_response(object)

    mss <- sapply(seq_len(nrow(object$submodel)), function (i) {
        if (is.na(rss[i]))  return (NA)

        x <- model.matrix(object, size = p[i], best = best[i])
        f <- stats_fitted(x, y, o, w)

        stats_mss(m, w, f, idf)
    })

    ll <- stats_log_lik(m, w, rss)

    list(m = m, n = n, idf = idf, p = p, rdf = m - p, rss = rss,
         mss = mss, ll = ll, sigma2_full = sigma2_full)
}


stats_lmSelect <- function (object, ...) {
    m <- object$nobs
    n <- object$nvar

    lm_full <- lm(object)
    rss_full <- deviance(lm_full)
    sigma2_full <- stats_sigma2(rss_full, m - n)

    idf <- if (object$intercept) 1 else 0

    o <- object$offset
    w <- object$weights

    best <- object$submodel$BEST
    ic <- object$submodel$IC

    y <- model_response(object)

    mss <- sapply(best, function (i) {
        if (is.na(ic[i]))  return (NA)

        x <- model.matrix(object, best = best[i])
        f <- stats_fitted(x, y, o, w)

        stats_mss(m, w, f, idf)
    })

    p <- object$submodel$SIZE
    rss <- deviance(object, best = NULL, drop = TRUE, na.rm = FALSE)
    rss <- as.vector(rss)
    ll <- stats_log_lik(m, w, rss)

    list(m = m, n = n, idf = idf, p = p, rdf = m - p, rss = rss,
         mss = mss, ll = ll, sigma2_full = sigma2_full)
}


sigma_stats <- function (object, ...) {
    sqrt(stats_sigma2(object$rss, object$rdf))
}


aic_stats <- function (object, ...) {
    stats_aic(object$ll, 2, object$p + 1)
}


bic_stats <- function (object, ...) {
    stats_bic(object$ll, object$m, object$p + 1)
}


r2_stats <- function (object, ...) {
    stats_r2(object$rss, object$mss)
}


r2adj_stats <- function (object, ..., r2 = r2_stats(object, ...)) {
    stats_r2adj(object$m, r2, object$idf, object$rdf)
}


pval_stats <- function (object, ...) {
    stats_pval(object$p, object$rss, object$mss,
               object$idf, object$rdf)
}


cp_stats <- function (object, ...) {
    stats_cp(object$m, object$p, object$rss, object$sigma2_full)
}
