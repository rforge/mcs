library("lmSubsets")


expand <- function (...) {
    args <- rev(list(...))
    args$stringsAsFactors <- FALSE

    ans <- do.call(expand.grid, args)
    rev(ans)
}


## Initialize benchmark.
##
## Arguments:
##   name    - (character)
##   nrep    - (integer) number of repetitions
##   timeout - (integer) timeout in seconds
##   seed    - (integer)
##
## Result: (list)
##   benchmark
##
benchmark <- function (name = "bm", nrep = 1, timeout = Inf,
                       seed = NULL) {
    bm <- list()

    bm$name <- name
    bm$nrep <- nrep
    bm$timeout <- timeout
    bm$seed <- if (!is.null(seed)) seed else as.numeric(Sys.time())

    bm
}


## Save benchmark.
##
## Arguments:
##   bm   - (list) benchmark
##   file - (character) file name
##
save_benchmark <- function (bm, file = NULL) {
    if (is.null(file))  file <- paste0(bm$name, ".RData")

    save(bm, file = file)

    invisible(bm)
}


## Load benchmark.
##
## Arguments:
##   name - (character) benchmark name
##
## Result: (list)
##   benchmark
##
load_benchmark <- function (name = "bm") {
    file <- paste0(name, ".RData")

    if (!file.exists(file)) {
        return (NULL)
    }

    load(file = file)

    bm
}


## Distribution.
##
## Arguments:
##   MEAN    - (double[])
##   SD      - (double[])
##   RATE    - (double[])
##   MEANLOG - (double[])
##   SDLOG   - (double[])
##
## Result: (list)
##   benchmark
##
dist_benchmark <- function (bm, MEAN = NA, SD = NA, RATE = NA, MEANLOG = NA,
                            SDLOG = NA) {
    ## table DIST:  distribution parameters
    DIST <- expand(ID = NA, MEAN = MEAN, SD = SD, RATE = RATE,
                   MEANLOG = MEANLOG, SDLOG = SDLOG)
    DIST["ID"] <- seq.int(NROW(bm$DIST) + 1, length.out = NROW(DIST))
    bm$DIST <- rbind(bm$DIST, DIST)

    ## done
    bm
}



## Define base cases.
##
## Arguments:
##   bm        - (list) benchmark
##   NOBS      - (integer[]) number of observations
##   NVAR      - (integer[]) number of variables (regressors)
##   NTRUE     - (integer[]) number of regressors in "true" model
##   DIST_ID   - (integer[]) distribution ID
##   INTERCEPT - (logical[])
##
## Rval: (list)
##   benchmark
##
cases_benchmark <- function (bm, NOBS, NVAR, NTRUE, DIST_ID,
                             INTERCEPT = TRUE) {
    ## table BASE:  base cases
    BASE <- expand(ID = NA, NOBS = NOBS, NVAR = NVAR, NTRUE = NTRUE,
                   DIST_ID = DIST_ID, INTERCEPT = INTERCEPT)
    BASE["ID"] <- seq.int(NROW(bm$BASE) + 1, length.out = NROW(BASE))
    bm$BASE <- rbind(bm$BASE, BASE)

    ## table SIMU:  simulations
    SIMU <- expand(ID = NA, BASE_ID = BASE[, "ID"],
                   REP = seq_len(bm$nrep), SEED = NA, RSS_F = NA,
                   BIC_F = NA, WHICH_T = NA, RSS_T = NA, BIC_T = NA)
    SIMU["ID"] <- seq.int(NROW(bm$SIMU) + 1, length.out = NROW(SIMU))
    bm$SIMU <- rbind(bm$SIMU, SIMU)

    ## seeding
    set.seed(bm$seed)
    bm$SIMU["SEED"] <- sample.int(.Machine$integer.max, size = NROW(bm$SIMU))

    ## simulations
    for (simu_id in SIMU[, "ID"]) {
        simu <- simu_benchmark(bm, simu_id)

        simu_ix <- with(bm$SIMU, which(ID == simu_id))
        bm$SIMU[simu_ix, "RSS_F"] <- simu$rss_full
        bm$SIMU[simu_ix, "BIC_F"] <- simu$bic_full
        bm$SIMU[simu_ix, "WHICH_T"] <- paste0(as.integer(simu$true), collapse = "")
        bm$SIMU[simu_ix, "RSS_T"] <- simu$rss_true
        bm$SIMU[simu_ix, "BIC_T"] <- simu$bic_true
    }

    ## done
    bm
}


## Define 'lmSubsets' cases and runs.
##
## Arguments:
##   bm        - (list) benchmark
##   NMIN      - (integer) 'lmSubsets' arg
##   NMAX      - (integer) 'lmSubsets' arg
##   TOLERANCE - (numeric) 'lmSubsets' arg
##   cases     - (integer[]) relevant base cases
##
## Result: (list)
##   benchmark
##
lmSubsets_benchmark <- function (bm, NMIN = NA, NMAX = NA,
                                 TOLERANCE = NA, cases) {
    if (is.null(bm$BASE)) {
        stop ("missing BASE")
    }

    if (missing(cases)) {
        cases <- bm$BASE[, "ID"]
    }

    ## table LM_SUBSETS:  cases
    LM_SUBSETS <- expand(ID = NA, BASE_ID = cases, NMIN = NMIN,
                         NMAX = NMAX, TOLERANCE = TOLERANCE)
    LM_SUBSETS["ID"] <- seq.int(NROW(bm$LM_SUBSETS) + 1,
                                length.out = NROW(LM_SUBSETS))
    bm$LM_SUBSETS <- rbind(bm$LM_SUBSETS, LM_SUBSETS)

    ## table RUN:  runs
    RUN <- expand(ID = NA, WHAT = "LM_SUBSETS",
                  CASE_ID = LM_SUBSETS[, "ID"],
                  REP = seq_len(bm$nrep), EXECUTED = FALSE,
                  INTERRUPTED = NA, USER = NA, SYSTEM = NA,
                  ELAPSED = NA)
    RUN["ID"] <- seq.int(NROW(bm$RUN) + 1, length.out = NROW(RUN))
    bm$RUN <- rbind(bm$RUN, RUN)

    ## table VALUE:  values
    for (run_id in RUN[, "ID"]) {
        case_id <- with(RUN, CASE_ID[ID == run_id])
        base_id <- with(LM_SUBSETS, BASE_ID[ID == case_id])
        nvar <- with(bm$BASE, NVAR[ID == base_id])

        VALUE <- expand(ID = NA, RUN_ID = run_id,
                        RANK = seq_len(nvar), RSS = NA, BIC = NA,
                        WHICH = NA, HIT = NA)
        VALUE["ID"] <- seq.int(NROW(bm$VALUE) + 1, length.out = NROW(VALUE))
        bm$VALUE <- rbind(bm$VALUE, VALUE)
    }

    ## done
    bm
}


## Define 'lmSelect' cases and runs.
##
## Arguments:
##   bm        - (list) benchmark
##   TOLERANCE - (numeric) 'lmSelect' arg
##   cases     - (integer[]) relevant base cases
##
## Result: (list)
##   benchmark
##
lmSelect_benchmark <- function (bm, TOLERANCE = NA, cases) {
    if (is.null(bm$BASE)) {
        stop ("missing BASE")
    }

    if (missing(cases)) {
        cases <- bm$BASE[, "ID"]
    }

    ## table LM_SELECT:  cases
    LM_SELECT <- expand(ID = NA, BASE_ID = cases,
                        TOLERANCE = TOLERANCE)
    LM_SELECT["ID"] <- seq.int(NROW(bm$LM_SELECT) + 1,
                               length.out = NROW(LM_SELECT))
    bm$LM_SELECT <- rbind(bm$LM_SELECT, LM_SELECT)

    ## table RUN:  runs
    RUN <- expand(ID = NA, WHAT = "LM_SELECT",
                  CASE_ID = LM_SELECT[, "ID"], REP = seq_len(bm$nrep),
                  EXECUTED = FALSE, INTERRUPTED = NA, USER = NA,
                  SYSTEM = NA, ELAPSED = NA)
    RUN["ID"] <- seq.int(NROW(bm$RUN) + 1, length.out = NROW(RUN))
    bm$RUN <- rbind(bm$RUN, RUN)

    ## table VALUE:  values
    VALUE <- expand(ID = NA, RUN_ID = RUN[, "ID"], RANK = NA,
                    RSS = NA, BIC = NA, WHICH = NA, HIT = NA)
    VALUE["ID"] <- seq.int(NROW(bm$VALUE) + 1,
                           length.out = NROW(VALUE))
    bm$VALUE <- rbind(bm$VALUE, VALUE)

    ## done
    bm
}


## Execute benchmark.
##
## Arguments:
##   bm    - (list) benchmark
##   cases - (integer[]) base cases to run
##   force - (logical) force execution
##   save  - (logical) save benchmark
##
## Result: (list)
##   benchmark
##
exec_benchmark <- function (bm, cases, force = FALSE, save = TRUE) {
    message ("Running benchmark...")

    bm <- exec_lmSubsets_benchmark(bm, cases, force = force, save = save)
    bm <- exec_lmSelect_benchmark(bm, cases, force = force, save = save)

    bm
}


## Execute 'lmSubsets' cases.
##
## Arguments:
##   bm    - (list) benchmark
##   cases - (integer[]) bases cases to run
##   force - (logical) force execution
##   save  - (logical) save benchmark
##
## Result: (list)
##   benchmark
##
exec_lmSubsets_benchmark <- function (bm, cases, force = FALSE,
                                      save = TRUE) {
    message ("Running LM_SUBSETS...")

    if (missing(cases)) {
        cases <- bm$LM_SUBSETS[, "ID"]
    } else {
        cases <- with(bm$LM_SUBSETS, ID[BASE_ID %in% cases])
    }

    runs <- with(bm$RUN, which((WHAT == "LM_SUBSETS") &
                               (CASE_ID %in% cases)))
    for (run_ix in runs) {
        case_id <- bm$RUN[run_ix, "CASE_ID"]
        base_id <- with(bm$LM_SUBSETS, BASE_ID[ID == case_id])
        rep <- bm$RUN[run_ix, "REP"]

        base_id <- with(bm$LM_SUBSETS, BASE_ID[ID == case_id])
        message ("  lmSubsets ", case_id, ":", rep, " (case ", base_id, ")...")

        if (!force && bm$RUN[run_ix, "EXECUTED"]) {
            message ("  Skipping.")

            next
        }

        case_ix <- with(bm$LM_SUBSETS, which(ID == case_id))
        nmin <- bm$LM_SUBSETS[case_ix, "NMIN"]
        nmax <- bm$LM_SUBSETS[case_ix, "NMAX"]
        tolerance <- bm$LM_SUBSETS[case_ix, "TOLERANCE"]

        cl <- call("lmSubsets")
        if (!is.na(nmin))  cl$nmin <- nmin
        if (!is.na(nmax))  cl$nmax <- nmax
        if (!is.na(tolerance))  cl$tolerance <- tolerance

        simu_id <- with(bm$SIMU, ID[(BASE_ID == base_id) &
                                    (REP == rep)])
        simu <- simu_benchmark(bm, simu_id, cl)

        bm$RUN[run_ix, "EXECUTED"] <- TRUE
        bm$RUN[run_ix, "INTERRUPTED"] <- simu$interrupted
        bm$RUN[run_ix, "USER"] <- simu$user
        bm$RUN[run_ix, "SYSTEM"] <- simu$system
        bm$RUN[run_ix, "ELAPSED"] <- simu$elapsed

        run_id <- bm$RUN[run_ix, "ID"]

        values <- with(bm$VALUE, which(RUN_ID == run_id))
        for (val_ix in values) {
            rank <- bm$VALUE[val_ix, "RANK"]
            size <- rank + simu$intercept

            if (is.na(simu$value$rank[1, size]))  next

            rss <- deviance(simu$value, size = size)

            which <- simu$value$which[, 1, size]
            if (simu$intercept)  which <- which[-1]

            bm$VALUE[val_ix, "RSS"] <- rss
            bm$VALUE[val_ix, "BIC"] <- bic(simu$value, size = size)
            bm$VALUE[val_ix, "WHICH"] <- paste(as.integer(which), collapse = "")
            bm$VALUE[val_ix, "HIT"] <- all(which == simu$true)
        }

        if (save)  save_benchmark(bm)

        message ("  ... done.  (", bm$RUN[run_ix, "ELAPSED"], "s)")
    }

    bm
}


## Execute 'lmSelect' cases.
##
## Arguments:
##   bm    - (list) benchmark
##   cases - (integer[]) cases to run
##   force - (logical) force execution
##   save  - (logical) save benchmark
##
## Return: (list)
##   benchmark
##
exec_lmSelect_benchmark <- function (bm, cases, force = FALSE, save = TRUE) {
    message ("Running LM_SELECT...")

    if (!missing(cases)) {
        cases <- with(bm$LM_SELECT, ID[BASE_ID %in% cases])
    } else {
        cases <- bm$LM_SELECT[, "ID"]
    }

    runs <- with(bm$RUN, which((WHAT == "LM_SELECT") &
                               (CASE_ID %in% cases)))
    for (run_ix in runs) {
        case_id <- bm$RUN[run_ix, "CASE_ID"]
        rep <- bm$RUN[run_ix, "REP"]

        base_id <- with(bm$LM_SELECT, BASE_ID[ID == case_id])
        message ("  lmSelect ", case_id, ":", rep, " (case ", base_id, ")...")

        if (!force && bm$RUN[run_ix, "EXECUTED"]) {
            message ("  Skipping.")

            next
        }

        case_ix <- with(bm$LM_SELECT, which(ID == case_id))
        tolerance <- bm$LM_SELECT[case_ix, "TOLERANCE"]

        cl <- call("lmSelect")
        if (!is.na(tolerance))  cl$tolerance <- tolerance

        base_id <- bm$LM_SELECT[case_ix, "BASE_ID"]
        simu_id <- with(bm$SIMU, ID[(BASE_ID == base_id) &
                                    (REP == rep)])
        simu <- simu_benchmark(bm, simu_id, cl)

        bm$RUN[run_ix, "EXECUTED"] <- TRUE
        bm$RUN[run_ix, "INTERRUPTED"] <- simu$interrupted
        bm$RUN[run_ix, "USER"] <- simu$user
        bm$RUN[run_ix, "SYSTEM"] <- simu$system
        bm$RUN[run_ix, "ELAPSED"] <- simu$elapsed

        rss <- deviance(simu$value)

        which <- simu$value$which[, 1]
        if (simu$intercept)  which <- which[-1]

        run_id <- bm$RUN[run_ix, "ID"]

        val_ix <- with(bm$VALUE, which(RUN_ID == run_id))
        bm$VALUE[val_ix, "RANK"] <- sum(which)
        bm$VALUE[val_ix, "RSS"] <- rss
        bm$VALUE[val_ix, "BIC"] <- bic(simu$value)
        bm$VALUE[val_ix, "WHICH"] <- paste0(as.integer(which), collapse = "")
        bm$VALUE[val_ix, "HIT"] <- all(which == simu$true)

        if (save)  save_benchmark(bm)

        message ("  ... done.  (", bm$RUN[run_ix, "ELAPSED"], "s)")
    }

    bm
}


## Simulate and run.
##
## Arguments:
##   bm      - (list) benchmark
##   simu_id - (integer) simulation ID
##   cl      - (call)
##   x, y    - (logical)
##
## Result: (list)
##
simu_benchmark <- function (bm, simu_id, cl, x = FALSE, y = FALSE) {
    ret_x <- x;  x <- NULL
    ret_y <- y;  y <- NULL

    rdist <- function (n, MEAN = 0, SD = 1, RATE = 1, MEANLOG = 0,
                       SDLOG = 1) {
        if (!missing(MEAN) || !missing(SD)) {
            return (rnorm(n, mean = MEAN, sd = SD))
        }

        if (!missing(RATE)) {
            return (rexp(n, rate = RATE))
        }

        if (!missing(SDLOG) || !missing(MEANLOG)) {
            return (rlnorm(n, meanlog = MEANLOG, sdlog = SDLOG))
        }

        stop ("undefined distribution")
    }

    simu_ix <- with(bm$SIMU, which(ID == simu_id))

    ## seeding
    set.seed(bm$SIMU[simu_ix, "SEED"])

    ## args
    base_id <- bm$SIMU[simu_ix, "BASE_ID"]
    base_ix <- with(bm$BASE, which(ID == base_id))

    nobs <- bm$BASE[base_ix, "NOBS"]
    nvar <- bm$BASE[base_ix, "NVAR"]
    ntrue <- bm$BASE[base_ix, "NTRUE"]
    intercept <- bm$BASE[base_ix, "INTERCEPT"]
    dist_id <- bm$BASE[base_ix, "DIST_ID"]


    ## coefficients
    coefs <- rep_len(0, length.out = nvar)
    coefs[sample.int(nvar, size = ntrue)] <- 1

    ## error
    dist_ix <- with(bm$DIST, which(ID == dist_id))
    dist <- as.list(bm$DIST[dist_ix, ])
    dist <- dist[(names(dist) != "ID") & !is.na(dist)]
    e <- do.call(rdist, c(n = nobs, dist))

    ## independent variables
    x <- rnorm(nobs * nvar)
    dim(x) <- c(nobs, nvar)

    ## dependent variable
    y <- cbind(intercept, x) %*% c(1, coefs) + e

    ## answer
    ans <- list()

    ans$nobs <- nobs
    ans$nvar <- nvar
    ans$intercept <- intercept
    ans$true <- as.logical(coefs)

    ## call
    if (!missing(cl)) {
        cl$formula <- x
        cl$y <- y
        cl$intercept <- intercept

        tryCatch({
            setTimeLimit(elapsed = bm$timeout)

            t <- system.time(z <- eval(cl))
            t <- summary(t)

            ans$value <- z
            ans$interrupted <- z$.interrupted
            ans$user <- t["user"]
            ans$system <- t["system"]
            ans$elapsed <- t["elapsed"]
        }, finally = {
            setTimeLimit(elapsed = NULL)
        })
    }
    
    ## stats
    x_full <- cbind(intercept, x)
    r_full <- qr.resid(qr(x_full), y)
    rss_full <- sum(r_full^2)
    ll_full <- lmSubsets:::stats_log_lik(nobs, NULL, rss_full)
    bic_full <- lmSubsets:::stats_bic(ll_full, nobs, nvar + intercept + 1)

    ans$rss_full <- rss_full
    ans$bic_full <- bic_full

    x_true <- cbind(intercept, x[, ans$true])
    r_true <- qr.resid(qr(x_true), y)
    rss_true <- sum(r_true^2)
    ll_true <- lmSubsets:::stats_log_lik(nobs, NULL, rss_true)
    bic_true <- lmSubsets:::stats_bic(ll_true, nobs, ntrue + intercept + 1)

    ans$rss_true <- rss_true
    ans$bic_true <- bic_true

    if (ret_x)  ans$x <- x
    if (ret_y)  ans$y <- y

    ## done
    ans
}


## Summarize benchmark.
##
## Arguments:
##   bm - (list) benchmark
##
## Result: (list)
##   benchmark
##
summary_benchmark <- function (bm) {
    by_what <- by(bm$RUN, bm$RUN[, "WHAT"], function (WHAT) {
        by_case <- by(WHAT, WHAT[, "CASE_ID"], function (CASE) {
            runs <- with(CASE, which((EXECUTED == TRUE) &
                                     (INTERRUPTED == FALSE)))

            nrun <- length(runs)

            row <- data.frame(ID = NA, WHAT = NA, CASE_ID = NA,
                              NRUN = nrun, NHIT = NA, TM_AVG = NA,
                              TM_MIN = NA, TM_MAX = NA)

            if (nrun == 0)  return (row)

            tm <- CASE[runs, "ELAPSED"]

            tm_avg <- mean(tm)
            tm_min <- min(tm)
            tm_max <- max(tm)

            nhit <- with(bm$VALUE, sum(HIT[RUN_ID %in% CASE[runs, "ID"]]))

            row[, "NHIT"] <- nhit
            row[, "TM_AVG"] <- tm_avg
            row[, "TM_MIN"] <- tm_min
            row[, "TM_MAX"] <- tm_max

            row
        })

        blk <- data.frame()
        for (case_id in names(by_case)) {
            row <- by_case[[case_id]]
            row[, "CASE_ID"] <- as.integer(case_id)
            blk <- rbind(blk, row)
        }

        blk
    })

    for (what in names(by_what)) {
        blk <- by_what[[what]]
        blk[, "WHAT"] <- what
        bm$SUMMARY <- rbind(bm$SUMMARY, blk)
    }
    bm$SUMMARY["ID"] <- seq_len(NROW(bm$SUMMARY))

    bm
}



## utility functions (Q&D)



get_lm_benchmark <- function (bm, base_id, rep, which) {
    simu_ix <- with(bm$SIMU, {
        base::which((BASE_ID == base_id) & (REP == rep))
    })
    simu_id <- bm$SIMU[simu_ix, "ID"]

    simu <- simu_benchmark(bm, simu_id, x = TRUE, y = TRUE)

    if (missing(which)) {
        ## return full model
        which <- TRUE  # all
    } else if (is.logical(which) && which) {
        ## return "true" model
        which <- simu$true
    } else if (is.character(which)) {
        ## return submodel
        which <- unlist(strsplit(which, split = ""))
        which <- as.logical(as.integer(which))
    }

    intercept <- with(bm$BASE, INTERCEPT[ID == base_id])
    if (intercept) {
        ans <- lm(simu$y ~ simu$x[, which] + 1)
    } else {
        ans <- lm(simu$y ~ simu$x[, which] + 0)
    }

    ans <- list(simu_ix = simu_ix, lm = ans)
    ans$rss <- deviance(ans$lm)
    ans$bic <- BIC(ans$lm)

    ans
}


get_lmFull_benchmark <- function (bm, base_id, rep) {
    get_lm_benchmark(bm, base_id, rep)
}


get_lmTrue_benchmark <- function (bm, base_id, rep) {
    get_lm_benchmark(bm, base_id, rep, TRUE)
}


get_lmSubsets_benchmark <- function (bm, run_id, rank) {
    run_ix <- with(bm$RUN, which(ID == run_id))
    case_id <- bm$RUN[run_ix, "CASE_ID"]
    base_id <- with(bm$LM_SUBSETS, BASE_ID[ID == case_id])
    rep <- bm$RUN[run_ix, "REP"]

    if (missing(rank)) {
        rank <- with(bm$BASE, NTRUE[ID == base_id])
    }

    val_ix <- with(bm$VALUE, which((RUN_ID == run_id) &
                                   (RANK == rank)))

    ans <- get_lm_benchmark(bm, base_id, rep, bm$VALUE[val_ix, "WHICH"])

    ans$run_ix <- run_ix
    ans$case_id <- case_id
    ans$base_id <- base_id
    ans$val_ix <- val_ix

    ans
}


get_lmSelect_benchmark <- function (bm, run_id) {
    run_ix <- with(bm$RUN, which(ID == run_id))
    case_id <- bm$RUN[run_ix, "CASE_ID"]
    base_id <- with(bm$LM_SELECT, BASE_ID[ID == case_id])
    rep <- bm$RUN[run_ix, "REP"]

    val_ix <- with(bm$VALUE, which(RUN_ID == run_id))

    ans <- get_lm_benchmark(bm, base_id, rep, bm$VALUE[val_ix, "WHICH"])

    ans$run_ix <- run_ix
    ans$case_id <- case_id
    ans$base_id <- base_id
    ans$val_ix <- val_ix

    ans
}
