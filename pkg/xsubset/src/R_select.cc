/**
 * @file R_select.cc
 */
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
#include <tuple>

#include "mcs/subset/lm.hh"
#include "mcs/subset/select.hh"
#include "mcs/subset/select1.hh"
#include "mcs/subset/subset_table.hh"
#include "mcs/subset/subset_table1.hh"
#include "mcs/subset/criteria.hh"


void
do_select(const mcs::subset::lm<double, int>& lm,
          const int mark,
          const std::vector<double>& tau,
          const int prad,
          const int nbest,
          double* rss,
          int* which,
          unsigned long& nodes)
{
  const int nvar = lm.nvar();

  mcs::subset::subset_table<double, int> t =
    mcs::subset::select(lm, mark, tau, prad, nbest, nodes);

  t.sort();

  for (int n = 1; n <= nvar; ++n)
    {
      for (int i = 0; i < nbest; ++i)
        {
          std::tuple<double, std::vector<int> > x = t.get(n, i);
          rss[i] = std::get<0>(x);
          for (int k = 0; k < n; k++)
            which[std::get<1>(x)[k]] = 1;
          which += nvar;
        }
      rss += nbest;
    }
}


template<template<typename V,
                  typename S>
         class Criterion>
void
do_select1(const mcs::subset::lm<double, int>& lm,
           const int mark,
           const Criterion<double, int>& crit,
           const double tau,
           const int prad,
           const int nbest,
           double* aic,
           int* which,
           double* rss,
           unsigned long& nodes)
{
  const int nvar = lm.nvar();

  mcs::subset::subset_table1<double, int, Criterion> t =
    mcs::subset::select1(lm, mark, crit, tau, prad, nbest, nodes);

  t.sort();

  for (int i = 0; i < nbest; ++i)
    {
      std::tuple<double, int, std::vector<int>, double> x = t.get(i);
      aic[i] = std::get<0>(x);
      for (int k = 0; k < std::get<1>(x); k++)
        which[std::get<2>(x)[k]] = 1;
      which += nvar;
      rss[i] = std::get<3>(x);
    }
}


extern "C"
void
R_select(const int* const nobs, const int* const nvar,
	 const double* const xy, const int* const mark,
	 const double* const penalty, const double* const tolerance,
	 const int* const pradius, const int* const nbest,
	 double* const rss, double* const aic, int* const which,
	 int* const nodes)
{
  mcs::subset::lm<double, int> lm(*nobs, *nvar, xy);

  unsigned long xnodes = 0;

  if (*penalty == 0)
    {
      std::vector<double> tau(*nvar + 1, 0);
      std::transform(tolerance, tolerance + *nvar, tau.begin() + 1,
                     std::bind2nd(std::plus<double>(), 1));

      do_select(lm, *mark, tau, *pradius, *nbest, rss, which, xnodes);
    }
  else
    {
      double tau = *tolerance + 1;

      do_select1(lm, *mark, mcs::subset::aic<double, int>(*penalty),
                 tau, *pradius, *nbest, aic, which, rss, xnodes);
    }

  *nodes = xnodes;
}
