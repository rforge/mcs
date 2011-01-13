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
          double* value,
          int* which,
          unsigned long& nodes)
{
  const int nvar = lm.nvar();

  mcs::subset::subset_table<double, int> t =
    mcs::subset::select(lm, mark, tau, prad, nbest, nodes);

  t.sort();

  for (int j = 1; j <= nvar; ++j)
    {
      for (int i = 0; i < nbest; ++i)
        {
          std::tuple<double, std::vector<int> > x = t.get(j, i);
          value[i] = std::get<0>(x);
          for (int k = 0; k < j; k++)
            which[std::get<1>(x)[k]] = 1;
          which += nvar;
        }
      value += nbest;
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
           double* value,
           int* which,
           unsigned long& nodes)
{
  const int nvar = lm.nvar();

  mcs::subset::subset_table1<double, int, Criterion> t =
    mcs::subset::select1(lm, mark, crit, tau, prad, nbest, nodes);

  t.sort();

  for (int i = 0; i < nbest; ++i)
    {
      std::tuple<double, int, std::vector<int> > x = t.get(i);
      value[i] = std::get<0>(x);
      for (int k = 0; k < std::get<1>(x); k++)
        which[std::get<2>(x)[k]] = 1;
      which += nvar;
    }
}


extern "C"
void
R_select(const int* const nobs, const int* const nvar,
	 const double* const xy, const int* const mark,
         const char** const criterion, const double* const tolerance,
	 const int* const pradius, const int* const nbest,
         double* const value, int* const which, int* const nodes)
{
  mcs::subset::lm<double, int> lm(*nobs, *nvar, xy);

  unsigned long xnodes = 0;

  if (std::strcmp(*criterion, "rss") == 0)
    {
      std::vector<double> tau(*nvar + 1, 0);
      std::transform(tolerance, tolerance + *nvar, tau.begin() + 1,
                     std::bind2nd(std::plus<double>(), 1));

      do_select(lm, *mark, tau, *pradius, *nbest, value, which, xnodes);
    }
  else if (std::strcmp(*criterion, "aic") == 0)
    {
      double tau = *tolerance + 1;

      do_select1(lm, *mark, mcs::subset::aic<double, int>(), tau, *pradius,
                 *nbest, value, which, xnodes);
    }
  else if (std::strcmp(*criterion, "bic") == 0)
    {
      double tau = *tolerance + 1;

      do_select1(lm, *mark, mcs::subset::aic<double, int>(std::log(*nobs)),
                 tau, *pradius, *nbest, value, which, xnodes);
    }
  else if (std::strcmp(*criterion, "cp") == 0)
    {
      double tau = *tolerance + 1;

      do_select1(lm, *mark, mcs::subset::cp<double, int>(), tau, *pradius,
                 *nbest, value, which, xnodes);
    }
  *nodes = xnodes;
}
