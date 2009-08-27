#include <algorithm>

#include "mcs/xsubset/xsubset.hh"


extern "C"
void
xsubset_R(int*    nobs, int* nvar, double* ay,
          int*    mark, int* prad, double* tau,
          double* rsel, int* isel,
          int*    nvis)
{
  using namespace MCS::xsubset;

  t_size arg_nobs = *nobs;
  t_size arg_nvar = *nvar;
  t_size arg_mark = *mark;
  t_size arg_prad = *prad;
  t_size arg_nvis;

  t_dmatrix   arg_ay = t_dmatrix(arg_nobs, arg_nvar, ay);
  t_dbl_array arg_tau(arg_nvar);
  t_dbl_array arg_rsel(arg_nvar);
  t_size_array arg_isel(arg_nvar * (arg_nvar - 1) / 2);

  std::copy(tau, tau + arg_nvar - 1,
            arg_tau.begin() + 1);

  select(arg_ay,           arg_mark,
         arg_prad,         arg_tau.begin(),
         arg_rsel.begin(), arg_isel.begin(),
         arg_nvis);

  std::copy(arg_rsel.begin() + 1,    arg_rsel.end(), rsel);
  std::copy(arg_isel.begin(), arg_isel.end(), isel);
  *nvis = arg_nvis;
}
