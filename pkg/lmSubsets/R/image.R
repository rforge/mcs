image.lmSubsets <- function(x, best = 1, size = NULL, which = NULL,
                            main = "Subset selection", xlab = "", ylab = NULL,
                            xaxs = "i", yaxs = "i", cex = 0.9,
                            col = gray.colors(2), hilite = "BIC", hilite.hue = 0,
                            uline = hilite, srt = 45, gap = "white", ...) {
    hilite.colors <- function (n, hue, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL) {
        ans <- seq.int(from = start^gamma, to = end^gamma, length.out = n)^(1/gamma)
        ans <- hsv(hue, ans, 1, alpha)
        rev(ans)
    }

    ## which
    arg.which <- which;  which <- NULL
    which <- x$which[, best, ]
    if (!is.null(arg.which))  which <- which[arg.which, , drop = FALSE]
    if (!is.null(size))  which <- which[, as.character(size), drop = FALSE]
    which <- which[ , apply(!is.na(which), 2L, all), drop = FALSE]
    
    M <- nrow(which)
    N <- ncol(which)

    ## heatmap
    heatmap <- structure(col[2L - which], dim = c(M, N))

    ## highlight
    if (!is.null(hilite) || !identical(hilite, FALSE)) {
        if (identical(hilite, TRUE))  hilite <- "BIC"    
        if (is.character(hilite)) {
            hilite <- summary(x, penalty = hilite)$summary$val[best, ]
            hilite <- as.numeric(names(hilite))[which.min(hilite)]
        }
        y_hilite <- match(hilite, as.numeric(colnames(which)))
        if (!all(is.na(y_hilite))) {
            hilite.hue <- rep(hilite.hue, length.out = length(y_hilite))
            hilite.hue <- hilite.hue[!is.na(y_hilite)]
            y_hilite <- y_hilite[!is.na(y_hilite)]
            for (i in seq_along(y_hilite)) {
                heatmap[, y_hilite[i]] <- hilite.colors(2, hilite.hue[i])[2L - which[, y_hilite[i], drop = FALSE]]
            }
        }
    } else {
        hilite <- NULL
    }

    ## underline
    rlab <- rownames(which)
    if (!is.null(uline) || !identical(uline, FALSE)) {
        if (identical(uline, TRUE))  uline <- "BIC"    
        if (is.character(uline)) {
            uline <- summary(x, penalty = uline)$summary$val[best, ]
            uline <- as.numeric(names(uline))[which.min(uline)]
        }
        y_uline <- match(uline, as.numeric(colnames(which)))
        if (!all(is.na(y_uline))) {
            y_uline <- which(apply(which[, y_uline, drop = FALSE], 1L, any))
            rlab <- sapply(seq_along(rlab), function (i) {
                if(i %in% y_uline)
                    parse(text = paste("underline(", rlab[i], ")", sep = ""))
                else
                    parse(text = rlab[i])
            })
        }
    } else {
        uline <- NULL
    }

    ## empty plot
    if (is.null(ylab))  ylab <- sprintf("Number of regressors (best = %s)", best)
    plot(0, 0, xlim = c(0, M) + 0.5, ylim = c(N, 0) + 0.5, type = "n",
         axes = FALSE, xaxs = xaxs, yaxs = yaxs, xlab = xlab, ylab = ylab,
         main = main, ...)
    
    ## paint heatmap
    graphics::rect(row(which) - 0.5, col(which) - 0.5, row(which) + 0.5, col(which) + 0.5,
                   col = heatmap, border = gap)

    ## x-axis labels (variable names)
    text(1L:M, par("usr")[3] + 0.02 * N, labels = rlab, srt = srt, 
         adj = c(1, 1), xpd = TRUE, cex = cex)  

    ## y-axis labels (subset sizes)
    clab0 <- as.numeric(colnames(which))
    clab <- pretty(clab0)
    clab[1L] <- pmax(clab[1L], clab0[1L])
    clab[length(clab)] <- pmin(clab[length(clab)], clab0[length(clab0)])
    clab <- sort(unique(c(clab, uline, hilite)))
    axis(2, at = match(clab, clab0), labels = clab)

    ## box
    box()
}
