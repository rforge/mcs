/**
 * @file select1.hh"
 */
#ifndef MCS_SUBSET_SELECT1_HH
#define MCS_SUBSET_SELECT1_HH


#include <vector>

#include "lm.hh"
#include "subset_table1.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size,
             template<typename V,
                      typename S>
             class Criterion>
    subset_table1<Value, Size, Criterion>
    select1(const lm<Value, Size>& x,
            Size mark,
            const std::vector<Value>& tau,
            Size prad,
            Size nbest,
            unsigned long& nodes);


    template<typename Value,
             typename Size,
             template<typename V,
                      typename S>
             class Criterion>
    subset_table1<Value, Size, Criterion>
    select1(const lm<Value, Size>& x,
            Size mark,
            const Criterion<Value, Size>& crit,
            const std::vector<Value>& tau,
            Size prad,
            Size nbest,
            unsigned long& nodes);


  }

}


#include "select1.cc"
#endif
