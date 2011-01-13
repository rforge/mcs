/**
 * @file select.hh"
 */
#ifndef MCS_SUBSET_SELECT_HH
#define MCS_SUBSET_SELECT_HH


#include <vector>

#include "lm.hh"
#include "subset_table.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    subset_table<Value, Size>
    select(const lm<Value, Size>& x,
           Size mark,
           const std::vector<Value>& tau,
           Size prad,
           Size nbest,
           unsigned long& nodes);


  }

}


#include "select.cc"
#endif
