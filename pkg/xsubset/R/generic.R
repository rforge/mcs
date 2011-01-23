
xsubset <- function (object, ...) UseMethod("xsubset")

refit <- function (object, ...) UseMethod("refit")

model.response <- function (object, ...) UseMethod("model.response")
model.response.default <- function (object, ...) stats::model.response(object)

model.weights <- function (object, ...) UseMethod("model.weights")
model.weights.default <- function (object, ...) stats::model.weights(object)

model.offset <- function (object, ...) UseMethod("model.offset")
model.offset.default <- function (object, ...) stats::model.offset(object)

environment <- function (object, ...) UseMethod("environment")
environment.default <- function (object, ...) base::environment(object)
