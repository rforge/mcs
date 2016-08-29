image.lmSubsets <- function(x, best = 1, size = NULL, which = NULL,
                            hilite = "BIC", uline = hilite,
                            main = "Subset selection", xlab = "", ylab = NULL,
                            col = "gray30", tint = 0.8, gaps = "white",
                            hilite.col = "red", hilite.tint = 0.8,
                            xaxs = "i", yaxs = "i", cex = 0.9,
                            srt = 45, ...) {
    tint.col <- function (col, factor = 0.8) {
        rgb <- col2rgb(col)
        rgb <- rgb + (255 - rgb) * factor
        rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
    }

    ## which
    which.arg <- which;  which <- NULL
    which <- x$which[, best, ]
    if (!is.null(which.arg))  which <- which[which.arg, , drop = FALSE]
    if (!is.null(size))  which <- which[, as.character(size), drop = FALSE]
    which <- which[ , apply(!is.na(which), 2L, all), drop = FALSE]
    
    M <- nrow(which)
    N <- ncol(which)

    ## heatmap
    tint.arg <- tint;  tint <- NULL
    tint <- col;  tint[2] <- tint.col(tint[1], tint.arg)
    heatmap <- structure(tint[2L - which], dim = c(M, N))

    ## highlight
    if (!is.null(hilite) || !identical(hilite, FALSE)) {
        if (identical(hilite, TRUE))  hilite <- "BIC"    
        if (is.character(hilite)) {
            hilite <- summary(x, penalty = hilite)$summary$val[best, ]
            hilite <- as.numeric(names(hilite))[which.min(hilite)]
        }
        y_hilite <- match(hilite, as.numeric(colnames(which)))
        if (!all(is.na(y_hilite))) {
            hilite.col <- rep(hilite.col, length.out = length(y_hilite))[!is.na(y_hilite)]
            hilite.tint <- rep(hilite.tint, length.out = length(y_hilite))[!is.na(y_hilite)]
            y_hilite <- y_hilite[!is.na(y_hilite)]
            for (i in seq_along(y_hilite)) {
                tint <- hilite.col[i];  tint[2] <- tint.col(tint[1], hilite.tint[i])
                heatmap[, y_hilite[i]] <- tint[2L - which[, y_hilite[i], drop = FALSE]]
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
    graphics::rect(row(heatmap) - 0.5, col(heatmap) - 0.5, row(heatmap) + 0.5, col(heatmap) + 0.5,
                   col = heatmap, border = gaps)

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
