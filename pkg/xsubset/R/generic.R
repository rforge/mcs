

xsubset <- function (object, ...) {
  UseMethod("xsubset")
}


refit <- function (object, ...) UseMethod("refit")


model.response <- function (object, ...) UseMethod("model.response")
model.response.default <- function (object, ...) {
  stats::model.response(object)
}
model.response.lm <- function (object, ...) {
  mf <- model.frame(object)
  mf[[1]]
}


model.weights <- function (object, ...) UseMethod("model.weights")
model.weights.default <- function (object, ...) {
  stats::model.weights(object)
}
model.weights.lm <- function (object, ...) {
  object$weights
}


model.offset <- function (object, ...) UseMethod("model.offset")
model.offset.default <- function (object, ...) {
  stats::model.offset(object)
}
model.offset.lm <- function (object, ...) {
  object$offset
}


environment <- function (object, ...) UseMethod("environment")
environment.default <- function (object, ...) {
  base::environment(object)
}
environment.lm <- function (object, ...) {
  environment(object$terms)
}
