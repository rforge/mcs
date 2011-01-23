

## workhorse method
##
## Args:
##   nobs        - (integer)
##   nvar        - (integer)
##   x           - (numeric)
##   y           - (numeric)
##   mark        - (integer)
##   aic.penalty - (numeric)
##   tolerance   - (numeric[])
##   pradius     - (integer)
##   nbest       - (integer)
##
## Rval: (list)
##   nobs        - (integer)
##   nvar        - (integer)
##   mark        - (integer)
##   aic.penalty - (numeric)
##   tolerance   - (numeric[])
##   nbest       - (integer)
##   value       - (numeric[,]|numeric)
##   which       - (integer[,,]|numeric[,])
##   .nodes      - (integer)
##
##
## Forwards call to C library.  Arguments are not processed.
##
xselect <- function (nobs, nvar, x, y, mark, aic.penalty,
                     tolerance, pradius, nbest)
{
  ## call C function
  C_args <- list(## in
                 nobs        = as.integer(nobs),
                 nvar        = as.integer(nvar),
                 xy          = as.double(cbind(x, y)),
                 mark        = as.integer(mark),
                 aic.penalty = as.double(aic.penalty),
                 tolerance   = as.double(tolerance),
                 pradius     = as.integer(pradius),
                 nbest       = as.integer(nbest),
                 ## out
                 value = as.double(rep(0, nvar * nbest)),
                 which = as.logical(rep(0, nvar * nbest * nvar)),
                 nodes = integer(1))
  C_rval <- do.call(".C", c(name = "R_select", C_args))

  if (aic.penalty == 0) {
    C_rval$value <- array(C_rval$value, dim = c(nbest, nvar))
    C_rval$which <- array(C_rval$which, dim = c(nvar, nbest, nvar))
  } else {
    C_rval$value <- C_rval$value
    C_rval$which <- array(C_rval$which, dim = c(nvar, nbest))
  }

  C_rval$xy <- NA

  C_rval
}
