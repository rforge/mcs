/**
 * @file subset_table.cc
 */
#ifndef MCS_SUBSET_SUBSET_TABLE_CC
#define MCS_SUBSET_SUBSET_TABLE_CC


#include <limits>
#include <vector>
#include <tuple>
#include <algorithm>

#include "../core/numeric/vector.hh"

#include "node.hh"
#include "subset_table.hh"


#define SUBSET_TABLE subset_table<Value, Size>


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    bool
    SUBSET_TABLE::entry_comp::operator ()(const std::tuple<Value, std::vector<Size> >& x1,
                                          const std::tuple<Value, std::vector<Size> >& x2)
    {
      return std::get<0>(x1) < std::get<0>(x2);
    }


    template<typename Value,
             typename Size>
    SUBSET_TABLE::subset_table(Size nvar,
                               Size nbest) :
      nvar_(nvar),
      nbest_(nbest),
      tab_(nvar_ + 1, std::vector<std::tuple<Value, std::vector<Size> > >(nbest_, std::make_tuple(std::numeric_limits<Value>::max(), std::vector<Size>(nvar_, 0))))
    {
    }

  
    template<typename Value,
             typename Size>
    Value
    SUBSET_TABLE::max_rss(const Size n) const
    {
      return std::get<0>(tab_[n][0]);
    }


    template<typename Value,
             typename Size>
    void
    SUBSET_TABLE::report(const node<Value, Size>& x)
    {
      report(x.nvar_, x.s_, x.mark_, x.rz_.col(x.nvar_));
    }


    template<typename Value,
             typename Size>
    template<template<typename V,
                      typename S>
             class Derived>
    void
    SUBSET_TABLE::report(Size n,
                         const std::vector<Size>& s,
                         const Size k,
                         const mcs::core::numeric::vector_base<Value, Size, Derived>& z)
    {
      Value rss = 0;
      //for (Size j = n + 1; (--j) - k; )
      do
        {
          rss += mcs::core::util::math::sqr(z(n));

          if (std::get<0>(tab_[n][0]) > rss)
            {
              std::pop_heap(tab_[n].begin(), tab_[n].end(), entry_comp());
              std::get<0>(tab_[n][nbest_ - 1]) = rss;
              std::copy_n(s.begin(), n, std::get<1>(tab_[n][nbest_ - 1]).begin());
              std::push_heap(tab_[n].begin(), tab_[n].end(), entry_comp());
            }
        }
      while ((--n) - k);
    }


    template<typename Value,
             typename Size>
    void
    SUBSET_TABLE::sort()
    {
      for (Size n = 0; n <= nvar_; ++n)
        std::sort_heap(tab_[n].begin(), tab_[n].end(), entry_comp());
    }


    template<typename Value,
             typename Size>
    std::tuple<Value, std::vector<Size> >
    SUBSET_TABLE::get(Size n,
                      Size i) const
    {
      return tab_[n][i];
    }


  }
  
}


#undef SUBSET_TABLE


#endif
