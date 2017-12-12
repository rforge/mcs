
lmSubsets <- function (formula, ...)
    UseMethod("lmSubsets")

lmSelect <- function (formula, ...)
    UseMethod("lmSelect")


refit <- function (object, ...)
    UseMethod("refit")


model_response <- function (data, ...)
    UseMethod("model_response")

model_response.default <- function (data, type = "any", ...)
    stats::model.response(data = data, type = type)
