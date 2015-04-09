

mcsSubset <- function (object, ...)
  UseMethod("mcsSubset")


refit <- function (object, ...)
  UseMethod("refit")


model.response <- function (data, ...)
  UseMethod("model.response")

model.response.default <- function(data, ...)
  stats::model.response(data)


model.weights <- function (x, ...)
  UseMethod("model.weights")

model.weights.default <- function (x, ...)
  stats::model.weights(x)


model.offset <- function (x, ...)
  UseMethod("model.offset")

model.offset.default <- function (x, ...)
  stats::model.offset(x)


environment <- function (fun, ...)
  UseMethod("environment")

environment.default <- function (fun, ...)
  base::environment(fun)
