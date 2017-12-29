## Copyright 2018  Marc Hofmann and Achim Zeileis
##
## This file is part of 'lmSubsets'.
##
## 'lmSubsets' is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## 'lmSubsets' is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with 'lmSubsets'.  If not, see <http://www.gnu.org/licenses/>.



ic <- function (penalty, ...) {
    if (is.function(penalty)) {
        f <- penalty

        s <- evalq(substitute(penalty), parent.frame())
        s <- deparse(s)

        format <- function (...) {
            s
        }

        eval <- function (obj, ..., drop = TRUE) {
            val <- deviance(obj, ..., drop = FALSE)
            val <- within(val, {
                df <- SIZE + 1
                IC <- mapply(f, SIZE, RSS)

                rm(RSS)
            })

            if (drop) {
                val <- with(val, {
                    structure(IC, df = df, SIZE = SIZE, BEST = BEST)
                })

                if ((length(val) == 1) && missing(drop)) {
                    val <- as.vector(val)
                }
            }

            val
        }

        penalty <- function (...) {
            f
        }
    } else if (tolower(penalty) == "aic") {
        format <- function (...) {
            "AIC"
        }

        eval <- function (obj, ...) {
            AIC(obj, ...)
        }

        penalty <- function (...) {
            2.0
        }
    } else if (tolower(penalty) == "bic") {
        format <- function (...) {
            "BIC"
        }

        eval <- function (obj, ...) {
            BIC(obj, ...)
        }

        penalty <- function (nobs, ...) {
            log(nobs)
        }
    } else if (is.numeric(penalty)) {
        k_ <- penalty

        format <- function (...) {
            paste0("AIC (k = ", format_default(k_, ...), ")")
        }

        eval <- function (obj, ..., k) {
            AIC(obj, ..., k = k_)
        }

        penalty <- function (...) {
            k_
        }
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


penalty_ic <- function (x, ...) {
    x$penalty(...)
}
