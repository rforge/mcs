

format_default <- function (x, na.encode = TRUE, ...) {
    ans <- format(x, na.encode = TRUE, ...)

    if (!na.encode) {
        ans <- ifelse(is.na(x), NA, ans)
    }

    ans
}


format_pval <- function (x, na.encode = TRUE, ...) {
    ans <- format.pval(x, ...)

    if (!na.encode) {
        ans <- ifelse(is.na(x), NA, ans)
    }

    ans
}


format_ordinal <- function (x) {
    ind <- c("th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "th", "th", "th", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th")
    paste0(x, ind[x %% 100 + 1])
}


format_which <- function (mask, seq_ind = "-", seq_sep = ",") {
    w <- which(mask)
    if (length(w) < 1) {
        return ("")
    }

    d <- c(diff(w), -1)
    z <- character(0)

    state <- 1
    pos <- 1
    while (state > 0) {
        switch (state, {
            ## state = 1
            if (d[pos] < 0) {
                z <- c(z, w[pos])
                state <- -1
            } else if (d[pos] == 1) {
                z <- c(z, w[pos], seq_ind)
                state <- 2
            } else {
                z <- c(z, w[pos], seq_sep)
                state <- 1
            }
        }, {
            ## state = 2
            if (d[pos] < 0) {
                z <- c(z, w[pos])
                state <- -1
            } else if (d[pos] == 1) {
                state <- 2
            } else {
                z <- c(z, w[pos], seq_sep)
                state <- 1
            }
        })

        pos <- pos + 1
    }

    paste(z, collapse = "")
}
