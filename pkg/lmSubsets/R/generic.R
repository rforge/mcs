
lmSubsets <- function (formula, ...)
    UseMethod("lmSubsets")

lmSelect <- function (formula, ...)
    UseMethod("lmSelect")


refit <- function (object, ...)
    UseMethod("refit")


is_NA <- function (x, ...)
    UseMethod("is_NA")

is_NA.default <- function (x, ...)
    base::is.na(x)


model_response <- function (data, ...)
    UseMethod("model_response")

model_response.default <- function (data, type = "any", ...)
    stats::model.response(data = data, type = type)


log_lik <- function (object, ...)
    UseMethod("log_lik")

aic <- function (object, ...)
    UseMethod("aic")

bic <- function (object, ...)
    UseMethod("bic")
