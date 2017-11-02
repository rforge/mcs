source("benchmark.R")


NAME <- "bm-03-limassol"


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
    ## DIST_ID | MEAN |   SD | RATE | MEANLOG | SDLOG | Distribution
    ## =============================================================
    ##       1 |  0.0 | 0.05 |    - |       - |     - |       normal
    ##       2 |  0.0 | 0.10 |    - |       - |     - |       normal
    ##       3 |  0.0 | 0.50 |    - |       - |     - |       normal
    ##       4 |  0.0 | 1.00 |    - |       - |     - |       normal
    ##       5 |  0.0 | 5.00 |    - |       - |     - |       normal
    ##
    bm <- dist_benchmark(bm, SD = c(0.05, 0.10, 0.50, 1.00, 5.00))


    ## Setup datasets.
    ##
    ## DATA_ID | NOBS | NVAR | NTRUE | DIST_ID | INTERCEPT
    ## ===================================================
    ##     1-5 | 1000 |   20 |    10 |     1-5 |       yes
    ## ---------------------------------------------------
    ##    6-10 | 1000 |   40 |    20 |     1-5 |       yes
    ## ---------------------------------------------------
    ##   11-15 | 1000 |   60 |    30 |     1-5 |       yes
    ## ---------------------------------------------------
    ##   16-20 | 1000 |   80 |    40 |     1-5 |       yes
    ##
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 20, NTRUE = 10, DIST_ID = 1:5)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 40, NTRUE = 20, DIST_ID = 1:5)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 60, NTRUE = 30, DIST_ID = 1:5)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 80, NTRUE = 40, DIST_ID = 1:5)


    ## Setup cases.
    ##
    ## CASE_ID | WHAT       | DATA_ID | TOLERANCE | NMIN | NMAX
    ## ========================================================
    ##    1-15 | LM_SUBSETS |    1-15 |       0.0 |    - |    -
    ## --------------------------------------------------------
    ##   16-30 | LM_SUBSETS |    1-15 |       0.1 |    - |    -
    ## --------------------------------------------------------
    ##   31-45 | LM_SELECT  |    1-15 |       0.0 |    - |    -
    ## --------------------------------------------------------
    ##   46-60 | LM_SELECT  |    1-15 |       0.1 |    - |    -
    ##
    bm <- lmSubsets_benchmark(bm, DATA_ID = 1:15, TOLERANCE = 0.0)
    bm <- lmSubsets_benchmark(bm, DATA_ID = 1:15, TOLERANCE = 0.1)
    bm <- lmSelect_benchmark(bm, DATA_ID = 1:15, TOLERANCE = 0.0)
    bm <- lmSelect_benchmark(bm, DATA_ID = 1:15, TOLERANCE = 0.1)


    ## Setup cases.
    ##
    ## CASE_ID | WHAT       | DATA_ID | TOLERANCE | NMIN | NMAX
    ## ========================================================
    ##   61-65 | LM_SUBSETS |   16-20 |       0.0 |   30 |   50
    ## --------------------------------------------------------
    ##   66-70 | LM_SUBSETS |   16-20 |       0.1 |   30 |   50
    ## --------------------------------------------------------
    ##   71-75 | LM_SELECT  |   16-20 |       0.0 |    - |    -
    ## --------------------------------------------------------
    ##   76-80 | LM_SELECT  |   16-20 |       0.1 |    - |    -
    ##
    bm <- lmSubsets_benchmark(bm, DATA_ID = 16:20, TOLERANCE = 0.0, NMIN = 30, NMAX = 50)
    bm <- lmSubsets_benchmark(bm, DATA_ID = 16:20, TOLERANCE = 0.1, NMIN = 30, NMAX = 50)
    bm <- lmSelect_benchmark(bm, DATA_ID = 16:20, TOLERANCE = 0.0)
    bm <- lmSelect_benchmark(bm, DATA_ID = 16:20, TOLERANCE = 0.1)
}


## Run benchmark.
##
bm <- run_benchmark(bm)


## Summarize benchmark.
##
bm <- summary_benchmark(bm)


## Save benchmark.
##
save_benchmark(bm)
