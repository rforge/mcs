source("benchmark.R")


NAME <- "bm-02"


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
    ## DIST_ID | MEAN | SD       | RATE | MEANLOG | SDLOG | Distribution
    ## =================================================================
    ##    1-10 |  0.0 | 1.0-10.0 |    - |       - |     - |       normal
    ##
    bm <- dist_benchmark(bm, SD = as.double(1:10))


    ## Setup datasets.
    ##
    ## DATA_ID | NOBS | NVAR | NTRUE | DIST_ID | INTERCEPT
    ## ===================================================
    ##    1-10 |   75 |   60 |    30 |    1-10 |       yes
    ## ---------------------------------------------------
    ##   11-20 |  100 |   80 |    40 |    1-10 |       yes
    ##
    bm <- data_benchmark(bm, NOBS = 75, NVAR = 60, NTRUE = 30, DIST_ID = 1:10)
    bm <- data_benchmark(bm, NOBS = 100, NVAR = 80, NTRUE = 40, DIST_ID = 1:10)


    ## Setup cases.
    ##
    ## CASE_ID | WHAT       | DATA_ID | TOLERANCE | NMIN | NMAX
    ## ========================================================
    ##    1-10 | LM_SUBSETS |    1-10 |         - |    - |    -
    ## --------------------------------------------------------
    ##   11-20 | LM_SELECT  |    1-10 |         - |    - |    -
    ##
    bm <- lmSubsets_benchmark(bm, DATA_ID = 1:10)
    bm <- lmSelect_benchmark(bm, DATA_ID = 1:10)


    ## Setup cases.
    ##
    ## CASE_ID | WHAT       | DATA_ID | TOLERANCE | NMIN | NMAX
    ## ========================================================
    ##   21-30 | LM_SUBSETS |   11-20 |         - |   30 |   50
    ## --------------------------------------------------------
    ##   31-40 | LM_SELECT  |   11-20 |         - |    - |    -
    ##
    bm <- lmSubsets_benchmark(bm, DATA_ID = 11:20, NMIN = 30, NMAX = 50)
    bm <- lmSelect_benchmark(bm, DATA_ID = 11:20)
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
