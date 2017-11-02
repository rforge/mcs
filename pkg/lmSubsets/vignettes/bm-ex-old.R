source("benchmark-old.R")


## benchmark name
NAME <- "bm-ex-old"


## cases (1-4) to execute
## e.g.  EXEC_CASES = 3
EXEC_CASES <- 1:4


## Load existing benchmark.
##
bm <- load_benchmark(NAME)


## Initialize benchmark.
##
if (is.null(bm)) {


    ## New benchmark.
    ##
    ## Repetitions per case:  5
    ##
    bm <- benchmark(name = NAME, nrep = 5, seed = 1)


    ## Setup distributions.
    ##
    ## DIST_ID | MEAN |  SD | RATE | MEANLOG | SDLOG | Distribution
    ## ============================================================
    ##       1 |  0.0 | 1.0 |    - |       - |     - | normal
    ## ------------------------------------------------------------
    ##       2 |    - |   - |  1.0 |       - |     - | exponential
    ##
    bm <- dist_benchmark(bm, SD = 1.0)
    bm <- dist_benchmark(bm, RATE = 1.0)


    ## Setup cases.
    ##
    ## CASE_ID | NVAR | NOBS | NTRUE | DIST_ID
    ## =======================================
    ##       1 |   20 |   25 |    10 |       1
    ##       2 |   20 |   25 |    10 |       2
    ## ---------------------------------------
    ##       3 |   40 |   50 |    20 |       1
    ##       4 |   40 |   50 |    20 |       2
    ##
    bm <- cases_benchmark(bm, NOBS = 25, NVAR = 20, NTRUE = 10, DIST_ID = 1:2)
    bm <- cases_benchmark(bm, NOBS = 50, NVAR = 40, NTRUE = 20, DIST_ID = 1:2)


    ## Setup runs.
    ##
    ## WHAT       | TOLERANCE | CASE_ID | #runs
    ## ========================================
    ## LM_SUBSETS |       0.0 |       1 |     5
    ## LM_SUBSETS |       0.0 |       2 |     5
    ## ----------------------------------------
    ## LM_SUBSETS |       0.1 |       3 |     5
    ## LM_SUBSETS |       0.1 |       4 |     5
    ## ----------------------------------------
    ## LM_SELECT  |       0.0 |       1 |     5
    ## LM_SELECT  |       0.0 |       2 |     5
    ## ----------------------------------------
    ## LM_SELECT  |       0.1 |       3 |     5
    ## LM_SELECT  |       0.1 |       4 |     5
    ## ========================================
    ## Total:                                40
    ##
    ## Note:  Number of repetitions per case: nrep = 5
    ##
    bm <- lmSubsets_benchmark(bm, TOLERANCE = 0.0, cases = 1:4)
    bm <- lmSubsets_benchmark(bm, TOLERANCE = 0.1, cases = 1:4)
    bm <- lmSelect_benchmark(bm, TOLERANCE = 0.0, cases = 1:4)
    bm <- lmSelect_benchmark(bm, TOLERANCE = 0.1, cases = 1:4)


}


## Execute benchmark.
##
bm <- exec_benchmark(bm, cases = EXEC_CASES)


## Summarize benchmark.
##
bm <- summary_benchmark(bm)


## Save benchmark.
##
save_benchmark(bm)
