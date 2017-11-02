## quick and dirty


validate_lmSelect_benchmark <- function (bm, tol = 1e-05) {
    T_CASE <- subset(bm$LM_SELECT[c("CASE_ID", "TOLERANCE")], {
        is.na(TOLERANCE) | TOLERANCE == 0.0
    })
    T_CASE <- merge(bm$CASE[c("ID", "DATA_ID")], T_CASE["CASE_ID"],
                    by.x = "ID", by.y = "CASE_ID")

    T_CASE0 <- subset(bm$LM_SUBSETS[c("CASE_ID", "TOLERANCE")], {
        is.na(TOLERANCE) | (TOLERANCE == 0.0)
    })
    T_CASE0 <- merge(bm$CASE[c("ID", "DATA_ID")], T_CASE0["CASE_ID"],
                     by.x = "ID", by.y = "CASE_ID")

    T_VALID <- merge(T_CASE, T_CASE0, by = "DATA_ID", suffix = c("", "0"))
    T_VALID["CASE_ID"] <- T_VALID["ID"];  T_VALID["ID"] <- NULL
    T_VALID["CASE0_ID"] <- T_VALID["ID0"];  T_VALID["ID0"] <- NULL

    T_VALID <- T_VALID[rep(seq_len(nrow(T_VALID)), each = bm$nrep), ]
    T_VALID["REP"] <- rep(seq_len(bm$nrep), length.out = nrow(T_VALID))

    T_VALID <- merge(T_VALID, bm$RUN[c("ID", "CASE_ID", "REP")],
                     by = c("CASE_ID", "REP"))
    T_VALID["RUN_ID"] <- T_VALID["ID"];  T_VALID["ID"] <- NULL

    T_VALID <- merge(T_VALID, bm$RUN[c("ID", "CASE_ID", "REP")],
                     by.x = c("CASE0_ID", "REP"), by.y = c("CASE_ID", "REP"))
    T_VALID["RUN0_ID"] <- T_VALID["ID"];  T_VALID["ID"] <- NULL

    T_VALID <- merge(T_VALID, bm$VALUE[c("RUN_ID", "RANK")], by = "RUN_ID")

    T_VALID["BIC"] <- NA

    for (ix in seq_len(nrow(T_VALID))) {
        data_id <- T_VALID[ix, "DATA_ID"]
        case_id <- T_VALID[ix, "CASE_ID"]
        case0_id <- T_VALID[ix, "CASE0_ID"]
        rep <- T_VALID[ix, "REP"]

        nvar <- with(bm$DATA, NVAR[ID == data_id])

        m <- get_value_benchmark(bm, case_id = case_id, rep = rep)
        valid <- TRUE
        for (rk in 1:nvar) {
            m0 <- get_value_benchmark(bm, case_id = case0_id, rep = rep,
                                      rank = rk)
            valid <- valid && (m$bic < m0$bic + tol)
        }

        T_VALID[ix, "BIC"] <- m$bic
        T_VALID[ix, "VALID"] <- valid
    }

    T_VALID <- T_VALID[c("CASE_ID", "DATA_ID", "CASE0_ID", "REP", "RUN_ID",
                         "RUN0_ID", "RANK", "BIC", "VALID")]

    invisible(T_VALID)
}


validate_tolerance_benchmark <- function (bm, tol = 1e-05) {
    ans1 <- validate_tolerance_lmSubsets_benchmark(bm, tol)
    ans2 <- validate_tolerance_lmSelect_benchmark(bm, tol)

    ans1 + ans2
}


validate_tolerance_lmSubsets_benchmark <- function (bm, tol = 1e-05) {
    T_CASE <- subset(bm$LM_SUBSETS[c("CASE_ID", "TOLERANCE")], TOLERANCE > 0.0)
    T_CASE <- merge(bm$CASE[c("ID", "DATA_ID")], T_CASE,
                    by.x = "ID", by.y = "CASE_ID")

    T_CASE0 <- subset(bm$LM_SUBSETS[c("CASE_ID", "TOLERANCE")], {
        is.na(TOLERANCE) | (TOLERANCE == 0.0)
    })
    T_CASE0 <- merge(bm$CASE[c("ID", "DATA_ID")], T_CASE0["CASE_ID"],
                     by.x = "ID", by.y = "CASE_ID")

    T_VALID <- merge(T_CASE, T_CASE0, by = "DATA_ID", suffix = c("", "0"))
    T_VALID["CASE_ID"] <- T_VALID["ID"];  T_VALID["ID"] <- NULL
    T_VALID["CASE0_ID"] <- T_VALID["ID0"];  T_VALID["ID0"] <- NULL

    T_VALID <- merge(T_VALID, bm$SIMU[c("DATA_ID", "REP", "RSS_FUL")],
                     by = "DATA_ID")

    T_VALID <- merge(T_VALID, bm$RUN[c("ID", "CASE_ID", "REP")],
                     by = c("CASE_ID", "REP"))
    T_VALID["RUN_ID"] <- T_VALID["ID"];  T_VALID["ID"] <- NULL

    T_VALID <- merge(T_VALID, bm$RUN[c("ID", "CASE_ID", "REP")],
                     by.x = c("CASE0_ID", "REP"), by.y = c("CASE_ID", "REP"))
    T_VALID["RUN0_ID"] <- T_VALID["ID"];  T_VALID["ID"] <- NULL

    ans <- 0

    for (ix in seq_len(NROW(T_VALID))) {
        case_id <- T_VALID[ix, "CASE_ID"]
        rep <- T_VALID[ix, "REP"]
        tolerance <- T_VALID[ix, "TOLERANCE"]
        case0_id <- T_VALID[ix, "CASE0_ID"]
        rss_full <- T_VALID[ix, "RSS_FUL"]
        run_id <- T_VALID[ix, "RUN_ID"]
        run0_id <- T_VALID[ix, "RUN0_ID"]

        rss <- with(bm$VALUE, RSS[RUN_ID == run_id])
        rss0 <- with(bm$VALUE, RSS[RUN_ID == run0_id])

        ok <- ((rss - rss_full) <= (1 + tolerance) * (rss0 - rss_full) + tol)
        ok <- ok[!is.na(ok)]
        ok <- all(ok)

        message ("Case ", case_id, "[", rep, "]: ", ok)

        ans <- ans + 1
    }

    c(total = NROW(T_VALID), ok = ans)
}


validate_tolerance_lmSelect_benchmark <- function (bm, tol = 1e-05) {
    T_CASE <- subset(bm$LM_SELECT[c("CASE_ID", "TOLERANCE")], TOLERANCE > 0.0)
    T_CASE <- merge(bm$CASE[c("ID", "DATA_ID")], T_CASE,
                    by.x = "ID", by.y = "CASE_ID")

    T_CASE0 <- subset(bm$LM_SELECT[c("CASE_ID", "TOLERANCE")], {
        is.na(TOLERANCE) | (TOLERANCE == 0.0)
    })
    T_CASE0 <- merge(bm$CASE[c("ID", "DATA_ID")], T_CASE0["CASE_ID"],
                     by.x = "ID", by.y = "CASE_ID")

    T_VALID <- merge(T_CASE, T_CASE0, by = "DATA_ID", suffix = c("", "0"))
    T_VALID["CASE_ID"] <- T_VALID["ID"];  T_VALID["ID"] <- NULL
    T_VALIDn["CASE0_ID"] <- T_VALID["ID0"];  T_VALID["ID0"] <- NULL

    T_VALID <- merge(T_VALID, bm$SIMU[c("DATA_ID", "REP", "BIC_INF")],
                     by = "DATA_ID")

    T_VALID <- merge(T_VALID, bm$RUN[c("ID", "CASE_ID", "REP")],
                     by = c("CASE_ID", "REP"))
    T_VALID["RUN_ID"] <- T_VALID["ID"];  T_VALID["ID"] <- NULL

    T_VALID <- merge(T_VALID, bm$RUN[c("ID", "CASE_ID", "REP")],
                     by.x = c("CASE0_ID", "REP"), by.y = c("CASE_ID", "REP"))
    T_VALID["RUN0_ID"] <- T_VALID["ID"];  T_VALID["ID"] <- NULL

    ans <- 0

    for (ix in seq_len(NROW(T_VALID))) {
        case_id <- T_VALID[ix, "CASE_ID"]
        rep <- T_VALID[ix, "REP"]
        tolerance <- T_VALID[ix, "TOLERANCE"]
        case0_id <- T_VALID[ix, "CASE0_ID"]
        bic_inf <- T_VALID[ix, "BIC_INF"]
        run_id <- T_VALID[ix, "RUN_ID"]
        run0_id <- T_VALID[ix, "RUN0_ID"]

        bic <- with(bm$VALUE, BIC[RUN_ID == run_id])
        bic0 <- with(bm$VALUE, BIC[RUN_ID == run0_id])

        ok <- ((bic - bic_inf) <= (1 + tolerance) * (bic0 - bic_inf) + tol)

        message ("Case ", case_id, "[", rep, "]: ", ok)

        ans <- ans + 1
    }

    c(total = NROW(T_VALID), ok = ans)
}
