## heatmap
##
## Args:
##   x           - (lmSubsets)
##   best        - (integer)
##   size        - (integer[])
##   which       - (integer[]|character[])
##   hilite      - (integer|character)
##   uline       - (integer|character)
##   col         - (character)
##   tint        - (numeric)
##   gaps        - (character)
##   hilite.col  - (character[])
##   hilite.tint - (numeric[])
##   ...         - forwarded to 'plot' function
##   (graphical parameters)
##
## Rval:  (numeric[,])
##
image.lmSubsets <- function(x, best = 1, size = NULL, which = NULL,
                            hilite = "BIC", uline = hilite,
                            main = "Subset selection", xlab = "", ylab = NULL,
                            col = gray.colors(1), tint = 0.8, gaps = "white",
                            hilite.col = "red", hilite.tint = tint,
                            xaxs = "i", yaxs = "i", cex = 0.9,
                            srt = 45, ...) {
    force(uline)
    
    ## util
    col.tint <- function (col, factor = 0.8) {
        rgb <- col2rgb(col)
        rgb <- rgb + (255 - rgb) * factor
        rgb(rgb[1L], rgb[2L], rgb[3L], maxColorValue = 255)
    }

    ## heatmap
    heatmap <- t(x$which[, best, , drop = TRUE])

    colnames(heatmap)[x$include] <- paste("+", colnames(heatmap)[x$include])
    colnames(heatmap)[x$exclude] <- paste("-", colnames(heatmap)[x$exclude])

    if (!is.null(which))  heatmap <- heatmap[, which, , drop = FALSE]
    if (!is.null(size))  heatmap <- heatmap[as.character(size), , drop = FALSE]
    heatmap <- heatmap[!apply(is.na(heatmap), 1L, all), , drop = FALSE]

    M <- nrow(heatmap)
    N <- ncol(heatmap)

    hot <- col;  cold <- col.tint(col, tint)
    heatmap.temp <- array(c(hot, cold)[2L - heatmap], dim = c(M, N))

    ## highlight
    if (!is.null(hilite) && !identical(hilite, FALSE)) {
        if (identical(hilite, TRUE))  hilite <- "BIC"    
        if (is.character(hilite)) {
            hilite <- summary(x, penalty = hilite)$summary$val[best, ]
            hilite <- as.numeric(names(hilite))[which.min(hilite)]
        }

        hilite.ix <- match(hilite, rownames(heatmap))

        for (i in seq_along(hilite.ix)) {
            hot <- hilite.col[i];  cold <- col.tint(hilite.col[i], hilite.tint)
            heatmap.temp[(row(heatmap) == hilite.ix[i]) &  heatmap] <- hot
            heatmap.temp[(row(heatmap) == hilite.ix[i]) & !heatmap] <- cold
        }
    } else {
        hilite <- NULL
    }

    ## underline
    x.labels <- colnames(heatmap)
    if (!is.null(uline) && !identical(uline, FALSE)) {
        if (identical(uline, TRUE))  uline <- "BIC"    
        if (is.character(uline)) {
            uline <- summary(x, penalty = uline)$summary$val[best, ]
            uline <- as.numeric(names(uline))[which.min(uline)]
        }

        uline.ix <- match(uline, rownames(heatmap))

        x.labels <- sapply(seq_along(x.labels), function (j) {
            if (any(heatmap[uline.ix, j])) {
                as.expression(substitute(underline(x), list(x = x.labels[j])))
            } else {
                as.expression(substitute(          x , list(x = x.labels[j])))
            }
        })
    } else {
        uline <- NULL
    }

    ## empty plot
    if (is.null(ylab))  ylab <- sprintf("Number of regressors (best = %s)", best)
    plot(0, 0, xlim = c(0, N) + 0.5, ylim = c(M, 0) + 0.5, type = "n",
         axes = FALSE, xaxs = xaxs, yaxs = yaxs, xlab = xlab, ylab = ylab,
         main = main, ...)
    
    ## paint heatmap
    graphics::rect(col(heatmap) - 0.5, row(heatmap) - 0.5, col(heatmap) + 0.5, row(heatmap) + 0.5,
                   col = heatmap.temp, border = gaps)

    ## x-axis labels (variable names)
    text(1L:N, par("usr")[3] + 0.02 * M, labels = x.labels, srt = srt, 
         adj = c(1, 1), xpd = TRUE, cex = cex)  

    ## y-axis labels (subset sizes)
    y.labels <- as.numeric(rownames(heatmap))
    y.pretty <- pretty(y.labels)
    y.pretty[1L]               <- max(y.pretty[1L]              , y.labels[1L]              )
    y.pretty[length(y.pretty)] <- min(y.pretty[length(y.pretty)], y.labels[length(y.labels)])
    y.pretty <- sort(unique(c(y.pretty, uline, hilite)))
    axis(2, at = match(y.pretty, y.labels), labels = y.pretty)

    ## box
    box()
}
