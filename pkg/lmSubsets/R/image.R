

image.lmSubsets <- function (x, main = NULL, sub = NULL, xlab = NULL,
                             ylab = NULL, size = NULL, best = 1,
                             which = NULL, col = gray.colors(2),
                             hilite, hilite_col = heat.colors,
                             hilite_lab = quote(lab), pad_size = 3,
                             pad_best = 1, pad_which = 3,
                             axis_pos = -4, axis_tck = -4,
                             axis_lab = -10, ..., axes = TRUE,
                             ann = par("ann")) {
    object <- x;  x <- NULL

    ## DATA

    ## heatmap
    if (is.null(size))  size <- object$sizes
    if (is.null(best))  best <- seq_len(object$nbest)
    if (is.null(which))  which <- seq_len(object$nvar)

    heatmap <- object$which[which, best, size, drop = FALSE]

    .cnames <- dimnames(heatmap)[[1]]
    .cnames[object$include[which]] <- paste0("+", .cnames[object$include[which]])
    .cnames[object$exclude[which]] <- paste0("-", .cnames[object$exclude[which]])

    .rnames <- rep(dimnames(heatmap)[[3]], each = dim(heatmap)[2])

    dim(heatmap) <- c(dim(heatmap)[1], dim(heatmap)[2] * dim(heatmap)[3])
    dimnames(heatmap) <- list(.cnames, .rnames)

    heatmap <- t(heatmap)
    heatmap <- heatmap[!apply(is.na(heatmap), 1L, all), , drop = FALSE]

    ## highlight
    if (missing(hilite)) {
        hilite <- NULL
    } else if (is.null(hilite)) {
        hilite <- seq_len(nrow(heatmap))
    } else if (is.matrix(hilite)) {
        hilite[, 1] <- match(hilite[, 1], size)
        hilite[, 2] <- match(hilite[, 2], best)

        hilite <- (hilite[, 1] - 1) * length(best) + hilite[, 2]
        hilite <- na.omit(hilite)
    }

    ## PLOT

    ## new plot
    plot.new()

    ## pars
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))

    ## margins
    par(mar = c(5, 4, 4, 4) + 0.1)

    ## plot window
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

    ## color map
    col <- matrix(col, ncol = 2)
    col <- col[rep_len(seq_len(nrow(col)), nrow(heatmap)), , drop = FALSE]

    if (!is.null(hilite) && !is.null(hilite_col)) {
        if (is.function(hilite_col))  {
            hilite_col <- hilite_col(length(hilite))
        } else {
            hilite_col <- rep_len(hilite_col, length(hilite))
        }

        col[hilite, 1] <- hilite_col
    }

    col <- ifelse(heatmap, col[row(heatmap), 1], col[row(heatmap), 2])

    ## padding
    pad_which <- pix2usr(x = pad_which, carry = 0)
    pad_size <- pix2usr(y = pad_size, carry = 0)
    pad_best <- pix2usr(y = pad_best, carry = 0)

    group <- match(rownames(heatmap), unique(rownames(heatmap)))

    padcnt_which <- seq_len(ncol(heatmap)) - 1
    padcnt_which <- rep(padcnt_which, each = nrow(heatmap))

    padcnt_best <- cumsum(c(0, group[-nrow(heatmap)] == group[-1]))
    padcnt_best <- rep.int(padcnt_best, ncol(heatmap))

    padcnt_size <- cumsum(c(0, group[-nrow(heatmap)] != group[-1]))
    padcnt_size <- rep.int(padcnt_size, ncol(heatmap))

    ## coords
    w <- (1 - (ncol(heatmap) - 1) * pad_which) / ncol(heatmap)
    h <- (1 - (max(padcnt_best) * pad_best) - (max(padcnt_size) * pad_size)) / nrow(heatmap)

    x <- (col(heatmap) - 1) * w + padcnt_which * pad_which
    y <- 1 - (row(heatmap) - 1) * h - padcnt_best * pad_best - padcnt_size * pad_size

    ## plot
    rect(x, y - h, x + w, y, col = col, border = NA)

    ## axes
    if (axes) {
        axis_pos <- rep_len(axis_pos, 2)
        axis_tck <- rep_len(axis_tck, 2)
        axis_lab <- rep_len(axis_lab, 2)

        ## y axis
        right <- pix2usr(x = axis_pos[2], carry = +1)
        tick <- pix2usr(x = axis_tck[2], carry = +1)
        left <- right + tick

        top <- y[which(c(TRUE, group[-nrow(heatmap)] != group[-1]))]
        bottom <- c(top[-1] + pad_size, 0)

        segments(x0 = right, y0 = top, x1 = right, y1 = bottom, lend = 1, xpd = TRUE)
        segments(x0 = right, y0 = top, x1 = left, y1 = top, lend = 1, xpd = TRUE)

        right <- pix2usr(x = axis_lab[2], carry = +1)
        text(x = right, y = top, labels = unique(rownames(heatmap)), adj = c(1, 1),
             cex = 0.9, xpd = TRUE)

        ## x axis
        left <- x[which(row(heatmap) == 1)]
        right <- c(left[-1] - pad_which, 1)

        top <- pix2usr(y = axis_pos[1], carry = +1)
        tick <- pix2usr(y = axis_tck[1], carry = +1)
        bottom <- top + tick

        segments(x0 = left, y0 = top, x1 = right, y1 = top, lend = 1, xpd = TRUE)
        segments(x0 = right, y0 = top, x1 = right, y1 = bottom, lend = 1, xpd = TRUE)

        labs <- colnames(heatmap)

        labs <- sapply(seq_along(labs), function (j) {
            env <- list(lab = labs[j])
            if (any(heatmap[hilite, j]))
                ans <- do.call(substitute, list(hilite_lab, env))
            else
                ans <- substitute(lab, env)
            as.expression(ans)
        })

        top <- pix2usr(y = axis_lab[1], carry = +1)
        text(x = right, y = top, labels = labs, adj = c(1, 1), srt = 45,
             cex = 0.9, xpd = TRUE)
    }

    ## annotations
    if (ann) {
        if (is.null(main))  main <- "All subsets"
        if (is.null(sub)) {
            sub <- paste0(format_ordinal(best), collapse = ", ")
            sub <- paste0("Best = {", sub, "}")
        }
        if (is.null(ylab))  ylab <- if (length(best) > 1) quote(Size %*% Best) else quote(Size)

        title(main = main, sub = sub, xlab = xlab, ylab = ylab)
    }

    ## done
    invisible(object)
}


image.lmSelect <- function (x, main = NULL, sub = NULL, xlab = NULL,
                            ylab = NULL, best = NULL, which = NULL,
                            col = gray.colors(2), hilite,
                            hilite_col = heat.colors,
                            hilite_lab = quote(lab), pad_best = 2,
                            pad_which = 2, axis_pos = -4,
                            axis_tck = -4, axis_lab = -10, ...,
                            axes = TRUE, ann = par("ann")) {
    object <- x;  x <- NULL

    ## DATA

    ## heatmap
    if (is.null(best))  best <- seq_len(object$nbest)
    if (is.null(which))  which <- seq_len(object$nvar)

    best <- rev(best)
    heatmap <- object$which[which, best, drop = FALSE]

    .size <- apply(heatmap, 2L, sum)
    dimnames(heatmap)[[2]] <-paste0(dimnames(heatmap)[[2]], " (", .size, ")")

    .cnames <- dimnames(heatmap)[[1]]
    .cnames[object$include[which]] <- paste0("+", .cnames[object$include[which]])
    .cnames[object$exclude[which]] <- paste0("-", .cnames[object$exclude[which]])
    dimnames(heatmap)[[1]] <- .cnames

    heatmap <- t(heatmap)

    ## highlight
    if (missing(hilite)) {
        hilite <- NULL
    } else if (is.null(hilite)) {
        hilite <- rev(seq_along(best))
    } else {
        hilite <- match(hilite, best)
    }

    ## PLOT

    ## new plot
    plot.new()

    ## pars
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))

    ## margins
    par(mar = c(5, 4, 4, 4) + 0.1)

    ## plot window
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

    ## color map
    col <- matrix(col, ncol = 2)
    col <- col[rev(rep_len(seq_len(nrow(col)), nrow(heatmap))), , drop = FALSE]

    if (!is.null(hilite) && !is.null(hilite_col)) {
        if (is.function(hilite_col)) {
            hilite_col <- hilite_col(length(hilite))
        } else {
            hilite_col <- rep_len(hilite_col, length(hilite))
        }

        col[hilite, 1] <- hilite_col
    }

    col <- ifelse(heatmap, col[row(heatmap), 1], col[row(heatmap), 2])

    ## padding
    pad_which <- pix2usr(x = pad_which, carry = 0)
    pad_best <- pix2usr(y = pad_best, carry = 0)

    padcnt_which <- seq_len(ncol(heatmap)) - 1
    padcnt_which <- rep(padcnt_which, each = nrow(heatmap))

    padcnt_best <- seq_len(nrow(heatmap)) - 1
    padcnt_best <- rep.int(padcnt_best, ncol(heatmap))

    ## coords
    w <- (1 - (ncol(heatmap) - 1) * pad_which) / ncol(heatmap)
    h <- (1 - (nrow(heatmap) - 1) * pad_best) / nrow(heatmap)

    x <- (col(heatmap) - 1) * w + padcnt_which * pad_which
    y <- 1 - (row(heatmap) - 1) * h - padcnt_best * pad_best

    ## plot
    rect(x, y - h, x + w, y, col = col, border = NA)

    ## axes
    if (axes) {
        axis_pos <- rep_len(axis_pos, 2)
        axis_tck <- rep_len(axis_tck, 2)
        axis_lab <- rep_len(axis_lab, 2)

        ## y axis
        right <- pix2usr(x = axis_pos[2], carry = +1)
        tick <- pix2usr(x = axis_tck[2], carry = +1)
        left <- right + tick

        top <- y[seq_len(nrow(heatmap))]
        bottom <- c(top[-1] + pad_best, 0)

        segments(x0 = right, y0 = top, x1 = right, y1 = bottom, lend = 1, xpd = TRUE)
        segments(x0 = right, y0 = top, x1 = left, y1 = top, lend = 1, xpd = TRUE)

        right <- pix2usr(x = axis_lab[2], carry = +1)
        text(x = right, y = top, labels = unique(rownames(heatmap)), adj = c(1, 1),
             cex = 0.9, xpd = TRUE)

        ## x axis
        left <- x[which(row(heatmap) == 1)]
        right <- c(left[-1] - pad_which, 1)

        top <- pix2usr(y = axis_pos[1], carry = +1)
        tick <- pix2usr(y = axis_tck[1], carry = +1)
        bottom <- top + tick

        segments(x0 = left, y0 = top, x1 = right, y1 = top, lend = 1, xpd = TRUE)
        segments(x0 = right, y0 = top, x1 = right, y1 = bottom, lend = 1, xpd = TRUE)

        labs <- colnames(heatmap)

        labs <- sapply(seq_along(labs), function (j) {
            env <- list(lab = labs[j])
            if (any(heatmap[hilite, j]))
                ans <- do.call(substitute, list(hilite_lab, env))
            else
                ans <- substitute(lab, env)
            as.expression(ans)
        })

        top <- pix2usr(y = axis_lab[1], carry = +1)
        text(x = right, y = top, labels = labs, adj = c(1, 1), srt = 45,
             cex = 0.9, xpd = TRUE)
    }

    ## annotations
    if (ann) {
        if (is.null(main))  main <- "Best subsets"

        title(main = main, sub = sub, xlab = xlab, ylab = ylab)
    }
}
