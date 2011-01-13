/**
 * @file subset_table.hh
 */
#ifndef MCS_SUBSET_SUBSET_TABLE_HH
#define MCS_SUBSET_SUBSET_TABLE_HH


#include <vector>
#include <tuple>

#include "node.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    class subset_table
    {


    private:

      struct entry_comp
      {

        bool
        operator ()(const std::tuple<Value, std::vector<Size> >& x1,
                    const std::tuple<Value, std::vector<Size> >& x2);

      };


    private:

      Size nvar_;

      Size nbest_;

      std::vector<std::vector<std::tuple<Value, std::vector<Size> > > > tab_;


    public:

      subset_table(Size nvar,
                   Size nbest);

      Value
      max_rss(Size n) const;

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

      std::tuple<Value, std::vector<Size> >
      get(Size n,
          Size i) const;

    };
    

  }

}


#include "subset_table.cc"
#endif
