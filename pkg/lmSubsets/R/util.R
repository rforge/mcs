parent.env <- function (env = parent.frame()) {
    base::parent.env(env)
}


pix2usr <- function (x, y, carry) {
    if (!missing(x)) {
        if (!missing(y))  warning ("unused argument 'y' will be discarded")

        ans <- grconvertX(x, "device", "user") - grconvertX(0, "device", "user")
        if (!missing(carry))  ans <- carry_sign(carry, x, ans)
    } else if (!missing(y)) {
        ans <- grconvertY(y, "device", "user") - grconvertY(0, "device", "user")
        if (!missing(carry))  ans <- carry_sign(carry, y, ans)
    } else {
        stop ("please specify 'x' or 'y'")
    }

    ans
}


carry_sign <- function (carry, a, b) {
    if (carry > 0) {
        sign(a) * abs(b)
    } else if (carry < 0) {
        -sign(a) * abs(b)
    } else {
        abs(b)
    }
}
