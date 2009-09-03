.TAU_MAX <- .Machine$double.xmax


.First.lib <- function(libname, pkgname) {
  library.dynam("mcsxsubset",
                package=pkgname,
                lib.loc=libname)
}


xsubset <- function (a, y, icpt = TRUE, prad = 1, tau = 0,
                     vmin = NULL, vmax = NULL) {
  ay <- cbind(a, y)
  if (is.data.frame(ay)) {
    ay <- as.matrix(ay)
  }
  nobs <- nrow(ay)
  nvar <- ncol(ay)
  if (length(tau) == 1) {
    tau <- rep(tau, nvar - 1)
  }
  if (is.null(vmin)) {
    vmin <- 1
  }
  if (is.null(vmax)) {
    vmax <- nvar - 1
  }
  mark <- 0
  xnames <- colnames(a)
  yname <- colnames(y)
  if (icpt) {
    ay <- cbind(1, ay)
    nvar <- nvar + 1
    tau <- c(0, tau)
    vmin <- vmin + 1
    vmax <- vmax + 1
    mark <- 1
  }
  if (vmin > 1) {
    tau[1:(vmin - 1)] <- .TAU_MAX
  }
  if (vmax < (nvar - 1)) {
    tau[(vmax + 1):(nvar - 1)] <- .TAU_MAX
  }
  args.in <- list(nobs = as.integer(nobs),
                  nvar = as.integer(nvar),
                  ay   = as.numeric(ay),
                  mark = as.integer(mark),
                  prad = as.integer(prad),
                  tau  = as.numeric(tau))
  args.out <- list(rsel = numeric(nvar - 1),
                   isel = integer(nvar * (nvar - 1) / 2),
                   nvis = integer(1))
  args <- c(name = "xsubset_R", args.in, args.out)
  t <- system.time(x <- do.call(".C", args))
  rsel <- list()
  isel <- list()
  for(i in vmin:vmax)
    {
      rsel[[i - icpt]] <- x$rsel[i]
      j <- i * (i - 1) / 2
      isel[[i - icpt]] <- x$isel[(j + 1 + icpt):(j + i)] + 1 - icpt
    }
  list(nobs = nobs,        nvar = nvar - icpt,
       icpt = icpt,        prad = prad,
       vmin = vmin - icpt, vmax = vmax - icpt,
       rsel = rsel,        isel = isel,
       nvis = x$nvis,      time = t[1],
       xnames = xnames,    yname = yname)
}


xsubset.model <- function (x, size) {
  list(size   = size,
       rss    = x$rsel[[size]],
       xnames = x$xnames[x$isel[[size]]],
       yname  = x$yname)
}


xsubset.print <- function (x) {
  cat <- function (...) base::cat(..., sep="")
  catln <- function (...) cat(..., "\n")
  catln("nobs = ", x$nobs)
  catln("nvar = ", x$nvar)
  catln("icpt = ", x$icpt)
  catln("vmin = ", x$vmin)
  catln("vmax = ", x$vmax)
  catln("yname = ", x$yname)
  catln()
  catln("size, RSS, xnames")
  catln("-----------------")
  for (i in x$vmin:x$vmax) {
    m <- xsubset.model(x, i)
    cat(m$size, ", ")
    cat(m$rss, ", ")
    cat("[", paste(sort(m$xnames), collapse=", "), "]")
    catln()
  }
}
