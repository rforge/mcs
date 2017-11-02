#ifndef R_LM_SUBSETS_CC
#define R_LM_SUBSETS_CC



#include <string>
#include <tuple>  // std::tie
#include <vector>



#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>



#include "r_interrupt.inc"



#include "mcs/core/matrix.hh"

#include "mcs/subset/aic.hh"
#include "mcs/subset/abba.hh"
#include "mcs/subset/bba.hh"
#include "mcs/subset/dca.hh"
#include "mcs/subset/hbba.hh"
#include "mcs/subset/subset.hh"
#include "mcs/subset/table.hh"

#include "mcs/util/misc.hh"



using matrix_cspan = mcs::core::matrix<const double&>;

using mcs::subset::subset_all;
using mcs::subset::abba_all;
using mcs::subset::hbba_all;
using mcs::subset::bba_all;
using mcs::subset::dca_all;

using mcs::subset::subset_best;
using mcs::subset::abba_best;
using mcs::subset::hbba_best;
using mcs::subset::bba_best;
using mcs::subset::dca_best;

using mcs::subset::aic;

using mcs::subset::table_all;
using mcs::subset::table_best;

using mcs::util::to_ordinal;



namespace {

const std::string algo_dflt = "DFLT";
const std::string algo_abba = "abba";
const std::string algo_bba = "bba";
const std::string algo_dca = "dca";
const std::string algo_hbba = "hbba";

}



extern "C"
SEXP
lmSubsets(
    SEXP r_algo,
    SEXP r_xy,
    SEXP r_mark,
    SEXP r_tau,
    SEXP r_nbest,
    SEXP r_prad
)
{
    int protect_cnt = 0;


    if (!Rf_isNull(r_algo) && !Rf_isString(r_algo))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'algo' must be a character string");
    }

    const auto algo = [&r_algo]() -> std::string {
        if (!Rf_isNull(r_algo))
            return CHAR(STRING_ELT(r_algo, 0));
        return algo_dflt;
    }();


    if (!Rf_isMatrix(r_xy))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'xy' must be a numeric matrix");
    }

    if (!Rf_isReal(r_xy))
    {
        r_xy = Rf_protect(Rf_coerceVector(r_xy, REALSXP));
        ++protect_cnt;
    }

    const int* const xy_dims =
        INTEGER(Rf_coerceVector(Rf_getAttrib(r_xy, R_DimSymbol), INTSXP));
    const int m = xy_dims[0];
    const int n = xy_dims[1] - 1;

    if (m <= n)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'xy' (%d x %d) must be a tall (or square) matrix", m, n + 1);
    }

    matrix_cspan ay_mat(m, n + 1, REAL(r_xy));

    if (!Rf_isNumeric(r_mark))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'mark' must be numeric");
    }

    const int mark = Rf_asInteger(r_mark);

    if (mark < 0)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'mark' [%d] must be a non-negative integer", mark);
    }


    if (!Rf_isNumeric(r_tau))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'tau' must be numeric vector");
    }

    if (LENGTH(r_tau) != n)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'tau' (%d) must be of length %d", LENGTH(r_tau), n);
    }

    if (!Rf_isReal(r_tau))
    {
        r_tau = Rf_protect(Rf_coerceVector(r_tau, REALSXP));
        ++protect_cnt;
    }

    const std::vector<double> tau(REAL(r_tau), REAL(r_tau) + n);


    if (!Rf_isNumeric(r_nbest))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'nbest' must be numeric");
    }

    const int nbest = Rf_asInteger(r_nbest);

    if (nbest < 1)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'nbest' [%d] must be positive integer", nbest);
    }


    if (!Rf_isNumeric(r_prad))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'prad' must be numeric");
    }

    const int prad = Rf_asInteger(r_prad);

    if (prad < 0)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'prad' [%d] must be a non-negative integer", prad);
    }


    r_interrupt_setup();


    table_all<double> tab;
    int node_cnt = -1;

    if (algo == algo_dflt)
    {
        tab = subset_all<double>(ay_mat, mark, tau, nbest, prad);
    }
    else if (algo == algo_abba)
    {
        std::tie(tab, node_cnt) =
            abba_all<double>(ay_mat, mark, tau, nbest, prad);
    }
    else if (algo == algo_hbba)
    {
        std::tie(tab, node_cnt) =
            hbba_all<double>(ay_mat, mark, tau, nbest, prad);
    }
    else if (algo == algo_bba)
    {
        std::tie(tab, node_cnt) =
            bba_all<double>(ay_mat, mark, nbest, prad);
    }
    else if (algo == algo_dca)
    {
        std::tie(tab, node_cnt) =
            dca_all<double>(ay_mat, mark, nbest, prad);
    }
    else
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'algo' [%s]: unexpected value", algo.c_str());
    }


    SEXP r_rss_dim = Rf_protect(Rf_allocVector(INTSXP, 2));
    ++protect_cnt;

    INTEGER(r_rss_dim)[0] = nbest;
    INTEGER(r_rss_dim)[1] = n;

    SEXP r_rss = Rf_protect(Rf_allocArray(REALSXP, r_rss_dim));
    ++protect_cnt;

    SEXP r_which_dim = Rf_protect(Rf_allocVector(INTSXP, 3));
    ++protect_cnt;

    INTEGER(r_which_dim)[0] = n;
    INTEGER(r_which_dim)[1] = nbest;
    INTEGER(r_which_dim)[2] = n;

    SEXP r_which = Rf_protect(Rf_allocArray(LGLSXP, r_which_dim));
    ++protect_cnt;

    SEXP r_rank_dim = Rf_protect(Rf_allocVector(INTSXP, 2));
    ++protect_cnt;

    INTEGER(r_rank_dim)[0] = nbest;
    INTEGER(r_rank_dim)[1] = n;

    SEXP r_rank = Rf_protect(Rf_allocArray(INTSXP, r_rank_dim));
    ++protect_cnt;

    for (int i = 0; i < nbest; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            const auto& res = tab[j][i];

            if (res)
            {
                REAL(r_rss)[j * nbest + i] = res.key();
                INTEGER(r_rank)[j * nbest + i] = res.size();

                for (int k = 0; k < n; ++k)
                {
                    LOGICAL(r_which)[(i + j * nbest) * n + k] = FALSE;
                }

                for (int k = 0; k < res.size(); ++k)
                {
                    LOGICAL(r_which)[(i + j * nbest) * n + res[k]] = TRUE;
                }
            }
            else
            {
                REAL(r_rss)[j * nbest + i] = NA_REAL;
                INTEGER(r_rank)[j * nbest + i] = NA_INTEGER;

                for (int k = 0; k < n; ++k)
                {
                    LOGICAL(r_which)[(i + j * nbest) * n + k] = NA_LOGICAL;
                }
            }
        }
    }

    SEXP r_ans_names = Rf_protect(Rf_allocVector(STRSXP, 5));
    ++protect_cnt;

    SET_STRING_ELT(r_ans_names, 0, Rf_mkChar(".nodes"));
    SET_STRING_ELT(r_ans_names, 1, Rf_mkChar(".interrupted"));
    SET_STRING_ELT(r_ans_names, 2, Rf_mkChar("rss"));
    SET_STRING_ELT(r_ans_names, 3, Rf_mkChar("which"));
    SET_STRING_ELT(r_ans_names, 4, Rf_mkChar("rank"));

    SEXP r_ans = Rf_protect(Rf_allocVector(VECSXP, 5));
    ++protect_cnt;

    Rf_namesgets(r_ans, r_ans_names);
    SET_VECTOR_ELT(r_ans, 0, Rf_ScalarInteger(node_cnt));
    SET_VECTOR_ELT(r_ans, 1, Rf_ScalarLogical(r_interrupt_flag()));
    SET_VECTOR_ELT(r_ans, 2, r_rss);
    SET_VECTOR_ELT(r_ans, 3, r_which);
    SET_VECTOR_ELT(r_ans, 4, r_rank);


    Rf_unprotect(protect_cnt);

    return r_ans;
}



extern "C"
SEXP
lmSelect(
    SEXP r_algo,
    SEXP r_xy,
    SEXP r_mark,
    SEXP r_penalty_arg,
    SEXP r_tau,
    SEXP r_nbest,
    SEXP r_prad
)
{
    int protect_cnt = 0;


    if (!Rf_isNull(r_algo) && !Rf_isString(r_algo))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'algo' must be a character string");
    }

    const auto algo = [&r_algo]() -> std::string {
        if (!Rf_isNull(r_algo))
            return CHAR(STRING_ELT(r_algo, 0));
        return algo_dflt;
    }();


    if (!Rf_isMatrix(r_xy))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'xy' must be a numeric matrix");
    }

    if (!Rf_isReal(r_xy))
    {
        r_xy = Rf_protect(Rf_coerceVector(r_xy, REALSXP));
        ++protect_cnt;
    }

    const int* const xy_dims =
        INTEGER(Rf_coerceVector(Rf_getAttrib(r_xy, R_DimSymbol), INTSXP));
    const int m = xy_dims[0];
    const int n = xy_dims[1] - 1;

    if (m <= n)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'xy' (%d x %d) must be a tall (or square) matrix", m, n + 1);
    }

    matrix_cspan ay_mat(m, n + 1, REAL(r_xy));


    if (!Rf_isNumeric(r_mark))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'mark' must be numeric");
    }

    const int mark = Rf_asInteger(r_mark);

    if (mark < 0)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'mark' [%d] must be a non-negative integer", mark);
    }


    if (!Rf_isNumeric(r_penalty_arg) && !Rf_isFunction(r_penalty_arg))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'penalty' must be numeric or a function");
    }

    const auto cost_aic = aic<double>(Rf_asReal(r_penalty_arg), m);

    SEXP r_size_arg = Rf_protect(Rf_allocVector(INTSXP, 1));
    ++protect_cnt;

    SEXP r_rss_arg = Rf_protect(Rf_allocVector(REALSXP, 1));
    ++protect_cnt;

    // SET_TAG(r_size_arg, Rf_install("size"));
    // SET_TAG(r_rss_arg, Rf_install("rss"));

    SEXP r_fcall = Rf_protect(Rf_lang3(r_penalty_arg, r_size_arg, r_rss_arg));
    ++protect_cnt;

    const auto cost_Rf = [&r_fcall, &r_size_arg, &r_rss_arg](
        const int size,
        const double rss
    ) -> double {
        INTEGER(r_size_arg)[0] = size;
        REAL(r_rss_arg)[0] = rss;
        return Rf_asReal(Rf_eval(r_fcall, R_GlobalEnv));
    };


    if (!Rf_isNumeric(r_tau))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'tau' must be numeric vector");
    }

    const double tau = Rf_asReal(r_tau);


    if (!Rf_isNumeric(r_nbest))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'nbest' must be numeric");
    }

    const int nbest = Rf_asInteger(r_nbest);

    if (nbest < 1)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'nbest' [%d] must be positive integer", nbest);
    }


    if (!Rf_isNumeric(r_prad))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'prad' must be numeric");
    }

    const int prad = Rf_asInteger(r_prad);

    if (prad < 0)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'prad' [%d] must be a non-negative integer", prad);
    }


    r_interrupt_setup();


    table_best<double> tab;
    int node_cnt = -1;

    if (algo == algo_dflt)
    {
        tab = Rf_isNumeric(r_penalty_arg)?
            subset_best<double, decltype(cost_aic)>(
                ay_mat, mark, cost_aic, tau, nbest, prad):
            subset_best<double, decltype(cost_Rf)>(
                ay_mat, mark, cost_Rf, tau, nbest, prad);
    }
    else if (algo == algo_abba)
    {
        std::tie(tab, node_cnt) = Rf_isNumeric(r_penalty_arg)?
            abba_best<double, decltype(cost_aic)>(
                ay_mat, mark, cost_aic, tau, nbest, prad):
            abba_best<double, decltype(cost_Rf)>(
                ay_mat, mark, cost_Rf, tau, nbest, prad);
    }
    else if (algo == algo_hbba)
    {
        std::tie(tab, node_cnt) = Rf_isNumeric(r_penalty_arg)?
            hbba_best<double, decltype(cost_aic)>(
                ay_mat, mark, cost_aic, tau, nbest, prad):
            hbba_best<double, decltype(cost_Rf)>(
                ay_mat, mark, cost_Rf, tau, nbest, prad);
    }
    else if (algo == algo_bba)
    {
        std::tie(tab, node_cnt) = Rf_isNumeric(r_penalty_arg)?
            bba_best<double, decltype(cost_aic)>(
                ay_mat, mark, cost_aic, nbest, prad):
            bba_best<double, decltype(cost_Rf)>(
                ay_mat, mark, cost_Rf, nbest, prad);
    }
    else if (algo == algo_dca)
    {
        std::tie(tab, node_cnt) = Rf_isNumeric(r_penalty_arg)?
            dca_best<double, decltype(cost_aic)>(
                ay_mat, mark, cost_aic, nbest, prad):
            dca_best<double, decltype(cost_Rf)>(
                ay_mat, mark, cost_Rf, nbest, prad);
    }
    else
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'algo' [%s]: unexpected value", algo.c_str());
    }


    SEXP r_penalty_dim = Rf_protect(Rf_allocVector(INTSXP, 1));
    ++protect_cnt;

    INTEGER(r_penalty_dim)[0] = nbest;

    SEXP r_penalty = Rf_protect(Rf_allocArray(REALSXP, r_penalty_dim));
    ++protect_cnt;

    SEXP r_which_dim = Rf_protect(Rf_allocVector(INTSXP, 2));
    ++protect_cnt;

    INTEGER(r_which_dim)[0] = n;
    INTEGER(r_which_dim)[1] = nbest;

    SEXP r_which = Rf_protect(Rf_allocArray(LGLSXP, r_which_dim));
    ++protect_cnt;

    SEXP r_rank_dim = Rf_protect(Rf_allocVector(INTSXP, 1));
    ++protect_cnt;

    INTEGER(r_rank_dim)[0] = nbest;

    SEXP r_rank = Rf_protect(Rf_allocArray(INTSXP, r_rank_dim));
    ++protect_cnt;

    for (int i = 0; i < nbest; ++i)
    {
        const auto& res = tab[i];

        if (res)
        {
            REAL(r_penalty)[i] = res.key();
            INTEGER(r_rank)[i] = res.size();

            for (int j = 0; j < n; ++j)
            {
                LOGICAL(r_which)[i * n + j] = FALSE;
            }

            for (int j = 0; j < res.size(); ++j)
            {
                LOGICAL(r_which)[i * n + res[j]] = TRUE;
            }
        }
        else
        {
            REAL(r_penalty)[i] = NA_REAL;
            INTEGER(r_rank)[i] = NA_INTEGER;

            for (int j = 0; j < n; ++j)
            {
                LOGICAL(r_which)[i * n + j] = NA_LOGICAL;
            }
        }
    }

    SEXP r_ans_names = Rf_protect(Rf_allocVector(STRSXP, 5));
    ++protect_cnt;

    SET_STRING_ELT(r_ans_names, 0, Rf_mkChar(".nodes"));
    SET_STRING_ELT(r_ans_names, 1, Rf_mkChar(".interrupted"));
    SET_STRING_ELT(r_ans_names, 2, Rf_mkChar("penalty"));
    SET_STRING_ELT(r_ans_names, 3, Rf_mkChar("which"));
    SET_STRING_ELT(r_ans_names, 4, Rf_mkChar("rank"));

    SEXP r_ans = Rf_protect(Rf_allocVector(VECSXP, 5));
    ++protect_cnt;

    Rf_namesgets(r_ans, r_ans_names);
    SET_VECTOR_ELT(r_ans, 0, Rf_ScalarInteger(node_cnt));
    SET_VECTOR_ELT(r_ans, 1, Rf_ScalarLogical(r_interrupt_flag()));
    SET_VECTOR_ELT(r_ans, 2, r_penalty);
    SET_VECTOR_ELT(r_ans, 3, r_which);
    SET_VECTOR_ELT(r_ans, 4, r_rank);


    Rf_unprotect(protect_cnt);

    return r_ans;
}



#endif
