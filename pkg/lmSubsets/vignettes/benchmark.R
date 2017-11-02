library("lmSubsets")


data_frame <- function (...) {
    data.frame(..., stringsAsFactors = FALSE)
}


expand_frame <- function (...) {
    args <- rev(list(...))
    args$stringsAsFactors <- FALSE

    ans <- do.call(expand.grid, args)
    rev(ans)
}


frame_id <- function (df_new, df) {
    ans <- seq_len(NROW(df_new))

    if (is.null(df))  return (ans)

    ans + df[NROW(df), "ID"]
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
##   name - (character) benchmark name
##
## Result: (list)
##   benchmark
##
save_benchmark <- function (bm, name = NULL) {
    if (is.null(name))  name <- bm$name

    file <- paste0(bm$name, ".RData")
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
load_benchmark <- function (name) {
    file <- paste0(name, ".RData")

    if (!file.exists(file)) {
        return (NULL)
    }

    load(file = file)

    bm
}


## Remove benchmark.
##
## Arguments:
##   bm   - (list) benchmark
##   name - (character) benchmark name
##
## Result: (logical)
##
remove_benchmark <- function (bm, name = NULL) {
    if (is.null(name))  name <- bm$name

    file <- paste0(name, ".RData")

    if (!file.exists(file)) {
        return (FALSE)
    }

    file.remove(file = file)
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
    T_DIST <- expand_frame(ID = NA, MEAN = MEAN, SD = SD, RATE = RATE,
                           MEANLOG = MEANLOG, SDLOG = SDLOG)
    T_DIST[, "ID"] <- frame_id(T_DIST, bm$DIST)
    bm$DIST <- rbind(bm$DIST, T_DIST)

    ## done
    bm
}



## Define datasets.
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
data_benchmark <- function (bm, NOBS, NVAR, NTRUE, DIST_ID,
                             INTERCEPT = TRUE) {
    ## table DATA: artificial data
    T_DATA <- expand_frame(ID = NA, NOBS = NOBS, NVAR = NVAR, NTRUE = NTRUE,
                           DIST_ID = DIST_ID, INTERCEPT = INTERCEPT)
    T_DATA[, "ID"] <- frame_id(T_DATA, bm$DATA)
    bm$DATA <- rbind(bm$DATA, T_DATA)

    ## table SIMU:  simulations
    T_SIMU <- expand_frame(ID = NA, DATA_ID = T_DATA[, "ID"],
                           REP = seq_len(bm$nrep), SEED = NA, RSS_FUL = NA,
                           BIC_FUL = NA, BIC_INF = NA, WHICH_TRU = NA,
                           RSS_TRU = NA, BIC_TRU = NA)
    T_SIMU[, "ID"] <- frame_id(T_SIMU, bm$SIMU)
    bm$SIMU <- rbind(bm$SIMU, T_SIMU)

    ## seeding
    set.seed(bm$seed)
    bm$SIMU[, "SEED"] <- sample.int(.Machine$integer.max, size = NROW(bm$SIMU))

    ## simulations
    for (simu_id in T_SIMU[, "ID"]) {
        simu <- simu_benchmark(bm, simu_id)

        simu_ix <- with(bm$SIMU, which(ID == simu_id))
        bm$SIMU[simu_ix, "RSS_FUL"] <- simu$rss_full
        bm$SIMU[simu_ix, "BIC_FUL"] <- simu$bic_full
        bm$SIMU[simu_ix, "BIC_INF"] <- simu$bic_inf
        bm$SIMU[simu_ix, "WHICH_TRU"] <- paste0(as.integer(simu$true), collapse = "")
        bm$SIMU[simu_ix, "RSS_TRU"] <- simu$rss_true
        bm$SIMU[simu_ix, "BIC_TRU"] <- simu$bic_true
    }

    ## done
    bm
}


## Define 'lmSubsets' cases and runs.
##
## Arguments:
##   bm        - (list) benchmark
##   DATA_ID   - (integer[]) relevant datasets
##   NMIN      - (integer) 'lmSubsets' arg
##   NMAX      - (integer) 'lmSubsets' arg
##   TOLERANCE - (numeric) 'lmSubsets' arg
##
## Result: (list)
##   benchmark
##
lmSubsets_benchmark <- function (bm, DATA_ID, NMIN = NA, NMAX = NA,
                                 TOLERANCE = NA) {
    if (missing(DATA_ID)) {
        stop ("missing argument 'DATA_ID'")
    }

    TMP <- expand_frame(DATA_ID = DATA_ID, NMIN = NMIN, NMAX = NMAX,
                        TOLERANCE = TOLERANCE)

    ## table CASE
    T_CASE <- data_frame(ID = NA, WHAT = "LM_SUBSETS",
                         DATA_ID = TMP[, "DATA_ID"])
    T_CASE[, "ID"] <- frame_id(T_CASE, bm$CASE)
    bm$CASE <- rbind(bm$CASE, T_CASE)

    ## table LM_SUBSETS
    T_LM_SUBSETS <- data_frame(ID = NA, CASE_ID = T_CASE[, "ID"],
                               NMIN = TMP[, "NMIN"], NMAX = TMP[, "NMAX"],
                               TOLERANCE = TMP[, "TOLERANCE"])
    T_LM_SUBSETS[, "ID"] <- frame_id(T_LM_SUBSETS, bm$LM_SUBSETS)
    bm$LM_SUBSETS <- rbind(bm$LM_SUBSETS, T_LM_SUBSETS)

    ## table RUN
    T_RUN <- expand_frame(ID = NA, CASE_ID = T_CASE[, "ID"],
                          REP = seq_len(bm$nrep), EXECUTED = FALSE,
                          INTERRUPTED = NA, USER = NA, SYSTEM = NA,
                          ELAPSED = NA)
    T_RUN[, "ID"] <- frame_id(T_RUN, bm$RUN)
    bm$RUN <- rbind(bm$RUN, T_RUN)

    ## table VALUE
    for (run_id in T_RUN[, "ID"]) {
        case_id <- with(T_RUN, CASE_ID[ID == run_id])
        data_id <- with(T_CASE, DATA_ID[ID == case_id])
        nvar <- with(bm$DATA, NVAR[ID == data_id])

        T_VALUE <- expand_frame(ID = NA, RUN_ID = run_id, RANK = seq_len(nvar),
                                RSS = NA, BIC = NA, WHICH = NA, HIT = NA)
        T_VALUE[, "ID"] <- frame_id(T_VALUE, bm$VALUE)
        bm$VALUE <- rbind(bm$VALUE, T_VALUE)
    }

    ## done
    bm
}


## Define 'lmSelect' cases and runs.
##
## Arguments:
##   bm        - (list) benchmark
##   DATA_ID   - (integer[]) relevant datasets
##   TOLERANCE - (numeric) 'lmSelect' arg
##
## Result: (list)
##   benchmark
##
lmSelect_benchmark <- function (bm, DATA_ID, TOLERANCE = NA) {
    if (missing(DATA_ID)) {
        stop ("missing argument 'DATA_ID'")
    }

    TMP <- expand_frame(DATA_ID = DATA_ID, TOLERANCE = TOLERANCE)

    ## table CASE
    T_CASE <- data_frame(ID = NA, WHAT = "LM_SELECT",
                         DATA_ID = TMP[, "DATA_ID"])
    T_CASE[, "ID"] <- frame_id(T_CASE, bm$CASE)
    bm$CASE <- rbind(bm$CASE, T_CASE)

    ## table LM_SELECT
    T_LM_SELECT <- data_frame(ID = NA, CASE_ID = T_CASE[, "ID"],
                              TOLERANCE = TMP[, "TOLERANCE"])
    T_LM_SELECT[, "ID"] <- frame_id(T_LM_SELECT, bm$LM_SELECT)
    bm$LM_SELECT <- rbind(bm$LM_SELECT, T_LM_SELECT)

    ## table RUN
    T_RUN <- expand_frame(ID = NA, CASE_ID = T_CASE[, "ID"],
                          REP = seq_len(bm$nrep), EXECUTED = FALSE,
                          INTERRUPTED = NA, USER = NA, SYSTEM = NA,
                          ELAPSED = NA)
    T_RUN[, "ID"] <- frame_id(T_RUN, bm$RUN)
    bm$RUN <- rbind(bm$RUN, T_RUN)

    ## table VALUE
    T_VALUE <- expand_frame(ID = NA, RUN_ID = T_RUN[, "ID"], RANK = NA,
                            RSS = NA, BIC = NA, WHICH = NA, HIT = NA)
    T_VALUE[, "ID"] <- frame_id(T_VALUE, bm$VALUE)
    bm$VALUE <- rbind(bm$VALUE, T_VALUE)

    ## done
    bm
}


## Run benchmark.
##
## Arguments:
##   bm    - (list) benchmark
##   what  - (character[]) restrict execution
##   cases - (integer[]) restrict execution
##   reps  - (integer[]) restrict execution
##   runs  - (integer[]) restrict execution
##   force - (logical) force execution
##   save  - (logical) save benchmark
##
## Result: (list)
##   benchmark
##
run_benchmark <- function (bm, what, cases, reps, runs, force = FALSE,
                           save = TRUE) {
    message ("Running benchmark '", bm$name, "'...")

    if (missing(what)) {
        what <- c("LM_SUBSETS", "LM_SELECT")
    }

    if (missing(cases)) {
        cases <- with(bm$CASE, ID[WHAT %in% what])
    }

    if (missing(reps)) {
        reps <- seq_len(bm$nrep)
    }

    if (missing(runs)) {
        runs <- with(bm$RUN, ID[(CASE_ID %in% cases) & (REP %in% reps)])
    }

    run_ixs <- with(bm$RUN, which(ID %in% runs))
    nrun <- NROW(run_ixs)

    for (run_ix in run_ixs) {
        run_id <- bm$RUN[run_ix, "ID"]
        case_id <- bm$RUN[run_ix, "CASE_ID"]
        rep <- bm$RUN[run_ix, "REP"]

        message ("  Run ", run_id, "/", nrun, ": case ", case_id, "[",
                 rep, "]...")

        if (!force && bm$RUN[run_ix, "EXECUTED"]) {
            message ("  ... skipping.")

            next
        }

        what <- with(bm$CASE, WHAT[ID == case_id])
        if (what == "LM_SUBSETS") {
            bm <- exec_lmSubsets_benchmark(bm, case_id, rep)
        } else {
            bm <- exec_lmSelect_benchmark(bm, case_id, rep)
        }

        if (save)  save_benchmark(bm)

        message ("  ... done.  (", bm$RUN[run_ix, "ELAPSED"], "s)")
    }

    message ("... done.  (benchmark '", bm$name, "')")

    bm
}


## Execute 'lmSubsets' case.
##
## Arguments:
##   bm      - (list) benchmark
##   case_id - (integer) case
##   rep     - (integer) repetition
##
## Result: (list)
##   benchmark
##
exec_lmSubsets_benchmark <- function (bm, case_id, rep = 1) {
    if (missing(case_id)) {
        stop ("missing argument 'case_id'")
    }

    case_ix <- with(bm$CASE, which(ID == case_id))
    data_id <- bm$CASE[case_ix, "DATA_ID"]

    lmSubsets_ix <- with(bm$LM_SUBSETS, which(CASE_ID == case_id))
    lmSubsets_id <- bm$LM_SUBSETS[lmSubsets_ix, "ID"]
    nmin <- bm$LM_SUBSETS[lmSubsets_ix, "NMIN"]
    nmax <- bm$LM_SUBSETS[lmSubsets_ix, "NMAX"]
    tolerance <- bm$LM_SUBSETS[lmSubsets_ix, "TOLERANCE"]

    cl <- call("lmSubsets")
    if (!is.na(nmin))  cl$nmin <- nmin
    if (!is.na(nmax))  cl$nmax <- nmax
    if (!is.na(tolerance))  cl$tolerance <- tolerance

    simu_id <- with(bm$SIMU, ID[(DATA_ID == data_id) & (REP == rep)])
    simu <- simu_benchmark(bm, simu_id, cl)

    run_ix <- with(bm$RUN, which((CASE_ID == case_id) & (REP == rep)))
    run_id <- bm$RUN[run_ix, "ID"]
    val_ixs <- with(bm$VALUE, which(RUN_ID == run_id))

    for (val_ix in val_ixs) {
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

    bm$RUN[run_ix, "EXECUTED"] <- TRUE
    bm$RUN[run_ix, "INTERRUPTED"] <- simu$interrupted
    bm$RUN[run_ix, "USER"] <- simu$user
    bm$RUN[run_ix, "SYSTEM"] <- simu$system
    bm$RUN[run_ix, "ELAPSED"] <- simu$elapsed

    bm
}


## Execute 'lmSelect' case.
##
## Arguments:
##   bm    - (list) benchmark
##   case_id - (integer) case
##   rep     - (integer) repetition
##
## Return: (list)
##   benchmark
##
exec_lmSelect_benchmark <- function (bm, case_id, rep = 1) {
    if (missing(case_id)) {
        stop ("missing argument 'case_id'")
    }

    case_ix <- with(bm$CASE, which(ID == case_id))
    data_id <- bm$CASE[case_ix, "DATA_ID"]

    lmSelect_ix <- with(bm$LM_SELECT, which(CASE_ID == case_id))
    lmSelect_id <- bm$LM_SELECT[lmSelect_ix, "ID"]
    tolerance <- bm$LM_SELECT[lmSelect_ix, "TOLERANCE"]

    cl <- call("lmSelect")
    if (!is.na(tolerance))  cl$tolerance <- tolerance

    simu_id <- with(bm$SIMU, ID[(DATA_ID == data_id) & (REP == rep)])
    simu <- simu_benchmark(bm, simu_id, cl)

    rss <- deviance(simu$value)

    which <- simu$value$which[, 1]
    if (simu$intercept)  which <- which[-1]

    run_ix <- with(bm$RUN, (CASE_ID == case_id) & (REP == rep))
    run_id <- bm$RUN[run_ix, "ID"]
    val_ix <- with(bm$VALUE, which(RUN_ID == run_id))

    bm$VALUE[val_ix, "RANK"] <- sum(which)
    bm$VALUE[val_ix, "RSS"] <- rss
    bm$VALUE[val_ix, "BIC"] <- bic(simu$value)
    bm$VALUE[val_ix, "WHICH"] <- paste0(as.integer(which), collapse = "")
    bm$VALUE[val_ix, "HIT"] <- all(which == simu$true)

    bm$RUN[run_ix, "EXECUTED"] <- TRUE
    bm$RUN[run_ix, "INTERRUPTED"] <- simu$interrupted
    bm$RUN[run_ix, "USER"] <- simu$user
    bm$RUN[run_ix, "SYSTEM"] <- simu$system
    bm$RUN[run_ix, "ELAPSED"] <- simu$elapsed

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
    data_id <- bm$SIMU[simu_ix, "DATA_ID"]
    data_ix <- with(bm$DATA, which(ID == data_id))

    nobs <- bm$DATA[data_ix, "NOBS"]
    nvar <- bm$DATA[data_ix, "NVAR"]
    ntrue <- bm$DATA[data_ix, "NTRUE"]
    intercept <- bm$DATA[data_ix, "INTERCEPT"]
    dist_id <- bm$DATA[data_ix, "DIST_ID"]


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

    ans$rss_full <- rss_full
    ans$rss_inf <- rss_full

    ans$bic_full <- lmSubsets:::stats_bic(ll_full, nobs, nvar + intercept + 1)
    ans$bic_inf <- lmSubsets:::stats_bic(ll_full, nobs, intercept + 1)

    x_true <- cbind(intercept, x[, ans$true])
    r_true <- qr.resid(qr(x_true), y)
    rss_true <- sum(r_true^2)
    ll_true <- lmSubsets:::stats_log_lik(nobs, NULL, rss_true)

    ans$rss_true <- rss_true
    ans$bic_true <- lmSubsets:::stats_bic(ll_true, nobs, ntrue + intercept + 1)

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
    bm$SUMMARY <- NULL

    for (case_id in bm$CASE[, "ID"]) {
        run_ixs <- with(bm$RUN, {
            which((CASE_ID %in% case_id) &
                  (EXECUTED == TRUE) &
                  (INTERRUPTED == FALSE))
        })

        nrun <- length(run_ixs)

        T_SUM1 <- data.frame(ID = NA, CASE_ID = case_id, NRUN = nrun,
                             NHIT = NA, TM_MIN = NA, TM_MAX = NA,
                             TM_AVG = NA)

        if (nrun > 0) {
            tm <- bm$RUN[run_ixs, "ELAPSED"]

            run_ids <- bm$RUN[run_ixs, "ID"]
            nhit <- with(bm$VALUE, sum(HIT[RUN_ID %in% run_ids]))

            T_SUM1[, "NHIT"] <- nhit
            T_SUM1[, "TM_AVG"] <- mean(tm)
            T_SUM1[, "TM_MIN"] <- min(tm)
            T_SUM1[, "TM_MAX"] <- max(tm)
        }

        bm$SUMMARY <- rbind(bm$SUMMARY, T_SUM1)
    }

    bm$SUMMARY[, "ID"] <- seq_len(NROW(bm$SUMMARY))

    bm
}



## extract models (quick and dirty)



get_simu_benchmark <- function (bm, data_id, rep, which) {
    simu_ix <- with(bm$SIMU, base::which((DATA_ID == data_id) & (REP == rep)))
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

    intercept <- with(bm$DATA, INTERCEPT[ID == data_id])
    if (intercept) {
        ans <- lm(simu$y ~ simu$x[, which] + 1)
    } else {
        ans <- lm(simu$y ~ simu$x[, which] + 0)
    }

    ans <- list(simu_id = simu_id, lm = ans)
    ans$rss <- deviance(ans$lm)
    ans$bic <- BIC(ans$lm)

    ans
}


get_full_benchmark <- function (bm, data_id, rep) {
    get_simu_benchmark(bm, data_id, rep)
}


get_true_benchmark <- function (bm, data_id, rep) {
    get_simu_benchmark(bm, data_id, rep, TRUE)
}


get_value_benchmark <- function (bm, case_id, rep, rank) {
    run_id <- with(bm$RUN, ID[(CASE_ID == case_id) & (REP == rep)])

    if (missing(rank)) {
        rank <- with(bm$VALUE, RANK[RUN_ID == run_id])
    }

    if (length(rank) > 1) {
        stop ("undetermined rank")
    }

    data_id <- with(bm$CASE, DATA_ID[ID == case_id])
    val_ix <- with(bm$VALUE, which((RUN_ID == run_id) & (RANK == rank)))

    ans <- get_simu_benchmark(bm, data_id, rep, bm$VALUE[val_ix, "WHICH"])

    ans$rep <- rep
    ans$rank <- rank
    ans$which <- bm$VALUE[val_ix, "WHICH"]
    ans$what <- with(bm$CASE, WHAT[ID == case_id])
    ans$run_id <- run_id
    ans$data_id <- data_id
    ans$val_id <- bm$VALUE[val_ix, "ID"]

    ans
}



## convert (quick and dirty)



convert_benchmark <- function (bm_old) {
    bm <- list()

    ## params
    bm$name <- bm_old$name
    bm$nrep <- bm_old$nrep
    bm$timeout <- bm_old$timeout
    bm$seed <- bm_old$seed

    ## distributions
    bm$DIST <- bm_old$DIST

    ## data
    bm$DATA <- bm_old$BASE

    ## simulations
    bm$SIMU <- expand_frame(ID = NA, DATA_ID = bm$DATA[, "ID"],
                            REP = seq_len(bm$nrep), SEED = NA, RSS_FUL = NA,
                            BIC_FUL = NA, BIC_INF = NA, WHICH_TRU = NA,
                            RSS_TRU = NA, BIC_TRU = NA)
    bm$SIMU[, "ID"] <- bm_old$SIMU[, "ID"]
    bm$SIMU[, "SEED"] <- bm_old$SIMU[, "SEED"]

    for (simu_id in bm$SIMU[, "ID"]) {

        simu <- simu_benchmark(bm, simu_id)

        simu_ix <- with(bm$SIMU, which(ID == simu_id))
        bm$SIMU[simu_ix, "RSS_FUL"] <- simu$rss_full
        bm$SIMU[simu_ix, "BIC_FUL"] <- simu$bic_full
        bm$SIMU[simu_ix, "BIC_INF"] <- simu$bic_inf
        bm$SIMU[simu_ix, "WHICH_TRU"] <- paste0(as.integer(simu$true), collapse = "")
        bm$SIMU[simu_ix, "RSS_TRU"] <- simu$rss_true
        bm$SIMU[simu_ix, "BIC_TRU"] <- simu$bic_true
    }

    ## cases
    nrun <- NROW(bm_old$RUN)

    T_CASE <- bm_old$RUN[seq.int(bm$nrep, nrun, bm$nrep),
                         c("ID", "WHAT", "CASE_ID")]
    T_CASE[, "ID"] <- T_CASE[, "ID"] / bm$nrep
    T_CASE["OLD_CASE_ID"] <- T_CASE["CASE_ID"];  T_CASE["CASE_ID"] <- NULL
    T_CASE["DATA_ID"] <- NA

    T_LM_SUBSETS <- bm_old$LM_SUBSETS;  T_LM_SUBSETS["CASE_ID"] <- NA
    T_LM_SELECT <- bm_old$LM_SELECT;  T_LM_SELECT["CASE_ID"] <- NA
    tmp <- list(LM_SUBSETS = T_LM_SUBSETS, LM_SELECT = T_LM_SELECT)

    for (ix in seq_len(NROW(T_CASE))) {
        what <- T_CASE[ix, "WHAT"]
        old_case_id <- T_CASE[ix, "OLD_CASE_ID"]
        old_what_ix <- with(bm_old[[what]], which(ID == old_case_id))
        old_base_id <- bm_old[[what]][old_what_ix, "BASE_ID"]

        T_CASE[ix, "DATA_ID"] <- old_base_id
        tmp[[what]][old_what_ix, "CASE_ID"] <- T_CASE[ix, "ID"]
    }

    bm$CASE <- T_CASE[, c("ID", "WHAT", "DATA_ID")]
    bm$LM_SUBSETS <- tmp$LM_SUBSETS[, c("ID", "CASE_ID", "NMIN", "NMAX",
                                        "TOLERANCE")]
    bm$LM_SELECT <- tmp$LM_SELECT[, c("ID", "CASE_ID", "TOLERANCE")]

    ## runs
    T_RUN <- bm_old$RUN
    T_RUN["OLD_CASE_ID"] <- T_RUN["CASE_ID"];  T_RUN["CASE_ID"] <- NULL

    T_RUN <- merge(T_RUN, T_CASE, by = c("WHAT", "OLD_CASE_ID"),
                   suffix = c("", ".case"))
    T_RUN["CASE_ID"] <- T_RUN["ID.case"]
    bm$RUN <- T_RUN[order(T_RUN[, "ID"]), {
        c("ID", "CASE_ID", "REP", "EXECUTED", "INTERRUPTED", "USER",
          "SYSTEM", "ELAPSED")
    }]

    ## values
    bm$VALUE <- bm_old$VALUE

    ## done
    bm
}
