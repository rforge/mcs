image.lmSubsets <- function(x, best = 1, size = NULL, which = NULL,
                            hilite = "BIC", uline = hilite,
                            main = "Subset selection", xlab = "", ylab = NULL,
                            col = "gray30", tint = 0.8, gaps = "white",
                            hilite.col = "red", hilite.tint = tint,
                            xaxs = "i", yaxs = "i", cex = 0.9,
                            srt = 45, ...) {
    ## args
    force(uline)
    force(hilite.tint)

    arg.which <- which;  which <- NULL
    arg.tint <- tint;  tint <- NULL

    ## util
    col.tint <- function (col, factor = 0.8) {
        rgb <- col2rgb(col)
        rgb <- rgb + (255 - rgb) * factor
        rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
    }

    ## which
    which <- x$which[, best, ]
    if (!is.null(arg.which))  which <- which[arg.which, , drop = FALSE]
    if (!is.null(size))  which <- which[, as.character(size), drop = FALSE]
    which <- which[ , apply(!is.na(which), 2L, all), drop = FALSE]
    
    M <- nrow(which)
    N <- ncol(which)

    ## heatmap
    tint <- col;  tint[2] <- col.tint(tint[1], arg.tint)
    heatmap <- structure(tint[2L - which], dim = c(M, N))

    ## highlight
    if (!is.null(hilite) || !identical(hilite, FALSE)) {
        if (identical(hilite, TRUE))  hilite <- "BIC"    
        if (is.character(hilite)) {
            hilite <- summary(x, penalty = hilite)$summary$val[best, ]
            hilite <- as.numeric(names(hilite))[which.min(hilite)]
        }
        hilite.ix <- match(hilite, as.numeric(colnames(which)))
        hilite.ok <- !is.na(hilite.ix)
        if (any(hilite.ok)) {
            hilite.col <- rep(hilite.col, length.out = length(hilite.ix))[hilite.ok]
            hilite.tint <- rep(hilite.tint, length.out = length(hilite.ix))[hilite.ok]
            hilite.ix <- hilite.ix[hilite.ok]
            for (i in seq_along(hilite.ix)) {
                tint <- hilite.col[i];  tint[2] <- col.tint(tint[1], hilite.tint[i])
                heatmap[, hilite.ix[i]] <- tint[2L - which[, hilite.ix[i], drop = FALSE]]
            }
        }
    } else {
        hilite <- NULL
    }

    ## underline
    x.names <- rownames(which)
    if (!is.null(uline) || !identical(uline, FALSE)) {
        if (identical(uline, TRUE))  uline <- "BIC"    
        if (is.character(uline)) {
            uline <- summary(x, penalty = uline)$summary$val[best, ]
            uline <- as.numeric(names(uline))[which.min(uline)]
        }
        uline.ix <- match(uline, as.numeric(colnames(which)))
        uline.ok <- !is.na(uline.ix)
        if (any(uline.ok)) {
            x.names <- ifelse(apply(which[, uline.ix, drop = FALSE], 1L, any),
                              paste("underline(", x.names, ")", sep = ""),
                              x.names)
            x.names <- sapply(x.names, function (x) parse(text = x))
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
    text(1L:M, par("usr")[3] + 0.02 * N, labels = x.names, srt = srt, 
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
