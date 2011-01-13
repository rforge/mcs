/**
 * @file select1.cc
 */
#ifndef MCS_SUBSET_SELECT1_CC
#define MCS_SUBSET_SELECT1_CC


#include <vector>

#include "lm.hh"
#include "subset_table1.hh"
#include "node_stack.hh"
#include "preorder.hh"
#include "node.hh"
#include "select1.hh"


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
            const Size mark,
            const std::vector<Value>& tau,
            const Size prad,
            const Size nbest,
            unsigned long& nodes)
    {
      Criterion<Value, Size> crit;

      return select1(x, mark, crit, tau, mark, prad, nbest, nodes);
    }



    template<typename Value,
             typename Size,
             template<typename V,
                      typename S>
             class Criterion>
    subset_table1<Value, Size, Criterion>
    select1(const lm<Value, Size>& x,
            const Size mark,
            const Criterion<Value, Size>& crit,
            const std::vector<Value>& tau,
            const Size prad,
            const Size nbest,
            unsigned long& nodes)
    {
      // criterion
      typename Criterion<Value, Size>::instance c = crit.for_model(x);

      // subset table
      subset_table1<Value, Size, Criterion> t(x.nvar(), nbest, c);

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
              if ((tau[n] * c.value(n, b)) >= t.max_val())
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
