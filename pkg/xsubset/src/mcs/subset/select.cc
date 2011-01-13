/**
 * @file select.cc
 */
#ifndef MCS_SUBSET_SELECT_CC
#define MCS_SUBSET_SELECT_CC


#include <vector>

#include "lm.hh"
#include "subset_table.hh"
#include "node_stack.hh"
#include "preorder.hh"
#include "node.hh"
#include "select.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    subset_table<Value, Size>
    select(const lm<Value, Size>& x,
           const Size mark,
           const std::vector<Value>& tau,
           const Size prad,
           const Size nbest,
           unsigned long& nodes)
    {
      // subset table
      subset_table<Value, Size> t(x.nvar(), nbest);

      // node stack
      node_stack<Value, Size> s(x.nvar());

      // preorder
      preorder<Value, Size> p(x.nvar(), mark, prad);
      
      // workspace (root)
      node<Value, Size> w(x, mark);

      // loop
      s.push(w);
      while (!s.empty())
        {
          // pop node
          s.pop(w);
   
          // preorder
          p.apply(w);
         
          // report
          t.report(w);
          
          // cut & spawn
          const Size k = w.mark();
          const Value b = w.bound();
          for (Size n = w.nvar(); (--n) - k; )
            {
              // cut
              if ((tau[n] * b) >= t.max_rss(n))
                continue;

              // spawn
              for (Size j = k; j < n; ++j)
                s.drop(w, j);
              break;
            }

          ++nodes;
        }

      // done
      return t;
    }


  }

}


#endif
