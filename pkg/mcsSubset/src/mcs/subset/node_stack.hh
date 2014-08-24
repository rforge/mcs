/**
 * @file node_stack.hh
 */
#ifndef MCS_SUBSET_NODE_STACK_HH
#define MCS_SUBSET_NODE_STACK_HH


#include <vector>

#include "node.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    class node_stack
    {

    private:

      std::vector<node<Value, Size> > stk_;

      typename std::vector<node<Value, Size> >::iterator bot_;

      typename std::vector<node<Value, Size> >::iterator top_;


    public:

      node_stack(Size nvar);

      void
      push(node<Value, Size>& x);

      void
      pop(node<Value, Size>& x);

      void
      drop(const node<Value, Size>& x,
           Size j);

      bool
      empty() const;

      Size
      size() const;

      Size
      capacity() const;

    };

  }

}


#include "node_stack.cc"
#endif
