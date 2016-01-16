
refit <- function (object, ...)
  UseMethod("refit")


model.response <- function (data, ...)
  UseMethod("model.response")

model.response.default <- function(data, type = "any", ...)
  stats::model.response(data = data, type = type)
