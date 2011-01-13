/**
 * @file subset_table1.cc
 */
#ifndef MCS_SUBSET_SUBSET1_TABLE_CC
#define MCS_SUBSET_SUBSET1_TABLE_CC


#include <limits>
#include <vector>
#include <tuple>
#include <algorithm>

#include "../core/numeric/vector.hh"

#include "node.hh"
#include "subset_table1.hh"


#define SUBSET_TABLE1 subset_table1<Value,       \
                                    Size,        \
                                    Criterion>


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size,
             template<typename V,
                      typename S>
             class Criterion>
    bool
    SUBSET_TABLE1::entry_comp::operator ()(const std::tuple<Value, Size, std::vector<Size> >& x1,
                                           const std::tuple<Value, Size, std::vector<Size> >& x2)
    {
      return std::get<0>(x1) < std::get<0>(x2);
    }


    template<typename Value,
             typename Size,
             template<typename V,
                      typename S>
             class Criterion>
    SUBSET_TABLE1::subset_table1(Size nvar,
                                 Size nbest,
                                 const typename Criterion<Value, Size>::instance& c) :
      nvar_(nvar),
      nbest_(nbest),
      c_(c),
      tab_(nbest_, std::make_tuple(std::numeric_limits<Value>::max(),
                                   Size(0), std::vector<Size>(nvar_, 0)))
    {
    }

  
    template<typename Value,
             typename Size,
             template<typename V,
                      typename S>
             class Criterion>
    Value
    SUBSET_TABLE1::max_val() const
    {
      return std::get<0>(tab_[0]);
    }


    template<typename Value,
             typename Size,
             template<typename V,
                      typename S>
             class Criterion>
    void
    SUBSET_TABLE1::report(const node<Value, Size>& x)
    {
      report(x.nvar_, x.s_, x.mark_, x.rz_.col(x.nvar_));
    }


    template<typename Value,
             typename Size,
             template<typename V,
                      typename S>
             class Criterion>
    template<template<typename V,
                      typename S>
             class Derived>
    void
    SUBSET_TABLE1::report(const Size n,
                          const std::vector<Size>& s,
                          const Size k,
                          const mcs::core::numeric::vector_base<Value, Size, Derived>& z)
    {
      Value rss = 0;
      for (Size j = n + 1; (--j) > k; )
        {
          rss += mcs::core::util::math::sqr(z(j));
          Value val = c_.value(j, rss);

          if (std::get<0>(tab_[0]) > val)
            {
              std::pop_heap(tab_.begin(), tab_.end(), entry_comp());
              std::get<0>(tab_[nbest_ - 1]) = val;
              std::get<1>(tab_[nbest_ - 1]) = j;
              std::copy_n(s.begin(), j, std::get<2>(tab_[nbest_ - 1]).begin());
              std::push_heap(tab_.begin(), tab_.end(), entry_comp());
            }
        }
    }


    template<typename Value,
             typename Size,
             template<typename V,
                      typename S>
             class Criterion>
    void
    SUBSET_TABLE1::sort()
    {
      std::sort_heap(tab_.begin(), tab_.end(), entry_comp());
    }


    template<typename Value,
             typename Size,
             template<typename V,
                      typename S>
             class Criterion>
    std::tuple<Value, Size, std::vector<Size> >
    SUBSET_TABLE1::get(Size i) const
    {
      return tab_[i];
    }


  }
  
}


#undef SUBSET_TABLE1


#endif
