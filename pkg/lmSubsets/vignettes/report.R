## quick and dirty


report_benchmark <- function (bm) {
    out <- NULL

    for (dist_id in bm$DIST[, "ID"]) {
        for (data_id in with(bm$DATA, ID[DIST_ID == dist_id])) {
            nvar <- with(bm$DATA, NVAR[ID == data_id])

            for (case_id in with(bm$CASE, ID[(WHAT == "LM_SUBSETS") &
                                             (DATA_ID == data_id)])) {
                nmin <- with(bm$LM_SUBSETS, NMIN[CASE_ID == case_id])
                nmax <- with(bm$LM_SUBSETS, NMAX[CASE_ID == case_id])
                tol <- with(bm$LM_SUBSETS, TOLERANCE[CASE_ID == case_id])

                tmin <- with(bm$SUMMARY, TM_MIN[CASE_ID == case_id])
                tmax <- with(bm$SUMMARY, TM_MAX[CASE_ID == case_id])
                tavg <- with(bm$SUMMARY, TM_AVG[CASE_ID == case_id])

                tmin <- lmSubsets:::format_numeric(tmin)
                tmax <- lmSubsets:::format_numeric(tmax)
                tavg <- lmSubsets:::format_numeric(tavg)

                out <- rbind(out, c(dist_id, "lmSubsets", data_id, nvar, nmin,
                                    nmax, tol, tmin, tmax, tavg))
            }
        }

        for (data_id in with(bm$DATA, ID[DIST_ID == dist_id])) {
            nvar <- with(bm$DATA, NVAR[ID == data_id])

            for (case_id in with(bm$CASE, ID[(WHAT == "LM_SELECT") &
                                             (DATA_ID == data_id)])) {
                nmin <- NA
                nmax <- NA
                tol <- with(bm$LM_SELECT, TOLERANCE[CASE_ID == case_id])

                tmin <- with(bm$SUMMARY, TM_MIN[CASE_ID == case_id])
                tmax <- with(bm$SUMMARY, TM_MAX[CASE_ID == case_id])
                tavg <- with(bm$SUMMARY, TM_AVG[CASE_ID == case_id])

                tmin <- lmSubsets:::format_numeric(tmin)
                tmax <- lmSubsets:::format_numeric(tmax)
                tavg <- lmSubsets:::format_numeric(tavg)

                out <- rbind(out, c(dist_id, "lmSelect", data_id, nvar, nmin,
                                    nmax, tol, tmin, tmax, tavg))
            }
        }
    }

    out <- ifelse(is.na(out), "", out)
    colnames(out) <- c("DIST_ID", "WHAT", "DATA_ID", "NVAR", "NMIN",
                       "NMAX", "TOL", "TMIN", "TMAX", "TAVG")
    rownames(out) <- rep("", nrow(out))
    print(out, quote = FALSE)

    invisible(bm)
}
