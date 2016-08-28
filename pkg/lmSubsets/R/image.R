image.lmSubsets <- function(x, best = 1, size = NULL, which = NULL,
  xlab = "", ylab = NULL, main = "Subset selection",
  xaxs = "i", yaxs = "i", cex = 0.9,
  col = gray.colors(2), rcol = 2, rect = "BIC", uline = rect,
  srt = 45, ...)
{
  ## all subset selections
  z <- x$which[, best, ]
  if(!is.null(which)) z <- z[which, , drop = FALSE]
  if(!is.null(size)) z <- z[, as.character(size), drop = FALSE]
  z <- z[ , apply(!is.na(z), 2L, all), drop = FALSE]

  ## empty plot
  nr <- nrow(z)
  nc <- ncol(z)
  if(is.null(ylab)) ylab <- sprintf("Number of regressors (best = %s)", best)
  plot(0, 0, xlim = c(0, nr) + 0.5, ylim = c(nc, 0) + 0.5, type = "n",
    axes = FALSE, xaxs = xaxs, yaxs = yaxs,
    xlab = xlab, ylab = ylab, main = main, ...)
  
  ## rectangles highlighting selection
  graphics::rect(row(z) - 0.5, col(z) - 0.5, row(z) + 0.5, col(z) + 0.5,
    col = col[2L - z], border = NA)

  ## highlight particular selection(s)
  if(!is.null(rect) | !identical(rect, FALSE)) {
    if(identical(rect, TRUE)) rect <- "BIC"    
    if(is.character(rect)) {
      rect <- summary(x, penalty = rect)$summary$val[best, ]
      rect <- as.numeric(names(rect))[which.min(rect)]
    }
    y_rect <- match(rect, as.numeric(colnames(z)))
    if(all(!is.na(y_rect))) {
      graphics::rect(rep(0.5, length(y_rect)), y_rect - 0.5,
        rep(nr, length(y_rect)) + 0.5, y_rect + 0.5, border = rcol)
    }
  } else {
    rect <- NULL
  }

  ## x-axis (variable) labels
  rlab <- rownames(z)
  if(!is.null(uline) | !identical(uline, FALSE)) {
    if(identical(uline, TRUE)) uline <- "BIC"    
    if(is.character(uline)) {
      uline <- summary(x, penalty = uline)$summary$val[best, ]
      uline <- as.numeric(names(uline))[which.min(uline)]
    }
    y_uline <- match(uline, as.numeric(colnames(z)))
    if(all(!is.na(y_uline))) {
      y_uline <- which(apply(z[, y_uline, drop = FALSE], 1L, any))
      rlab <- sapply(seq_along(rlab), function(i)
        if(i %in% y_uline) parse(text = paste("underline(", rlab[i], ")", sep = "")) else parse(text = rlab[i]))
    }
  } else {
    uline <- NULL
  }
  text(1L:nr, par("usr")[3] + 0.02 * nc, labels = rlab, srt = srt, 
    adj = c(1, 1), xpd = TRUE, cex = cex)  

  ## y-axis (subset) labels
  clab0 <- as.numeric(colnames(z))
  clab <- pretty(clab0)
  clab[1L] <- pmax(clab[1L], clab0[1L])
  clab[length(clab)] <- pmin(clab[length(clab)], clab0[length(clab0)])
  clab <- sort(unique(c(clab, uline, rect)))
  axis(2, at = match(clab, clab0), labels = clab)
  box()
}
