/**
 * @file subset_table1.hh
 */
#ifndef MCS_SUBSET_SUBSET1_TABLE_HH
#define MCS_SUBSET_SUBSET1_TABLE_HH


#include <vector>
#include <tuple>

#include "node.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size,
             template<typename V,
                      typename S>
             class Criterion>
    class subset_table1
    {


    private:

      struct entry_comp
      {

        bool
        operator ()(const std::tuple<Value, Size, std::vector<Size>, Value>& x1,
                    const std::tuple<Value, Size, std::vector<Size>, Value>& x2);

      };


    private:

      Size nvar_;

      Size nbest_;

      typename Criterion<Value, Size>::instance c_;

      std::vector<std::tuple<Value, Size, std::vector<Size>, Value> > tab_;


    public:

      subset_table1(Size nvar,
                    Size nbest,
                    const typename Criterion<Value, Size>::instance& c);

      Value
      max_val() const;

      void
      report(const node<Value, Size>& x);

      template<template<typename V,
                        typename S>
               class Derived>
      void
      report(Size n,
             const std::vector<Size>& s,
             Size k,
             const mcs::core::numeric::vector_base<Value, Size, Derived>& z);

      void
      sort();

      std::tuple<Value, Size, std::vector<Size>, Value>
      get(Size i) const;

    };
    

  }

}


#include "subset_table1.cc"
#endif
