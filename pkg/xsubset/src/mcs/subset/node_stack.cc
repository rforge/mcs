/**
 * @file node_stack.cc
 */
#ifndef MCS_SUBSET_NODE_STACK_CC
#define MCS_SUBSET_NODE_STACK_CC


#include "node.hh"
#include "node_stack.hh"


#define NODE_STACK node_stack<Value, Size>


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    NODE_STACK::node_stack(const Size nvar) :
      stk_(nvar + 1, node<Value, Size>(nvar)),
      bot_(stk_.begin()),
      top_(stk_.begin())
    {
    }


    template<typename Value,
             typename Size>
    void
    NODE_STACK::push(node<Value, Size>& x)
    {
      x.swap(*(top_++));
    }


    template<typename Value,
             typename Size>
    void
    NODE_STACK::pop(node<Value, Size>& x)
    {
      x.swap(*(--top_));
    }


    template<typename Value,
             typename Size>
    void
    NODE_STACK::drop(const node<Value, Size>& x,
                     const Size j)
    {
      x.drop(j, *(top_++));
    }


    template<typename Value,
             typename Size>
    bool
    NODE_STACK::empty() const
    {
      return bot_ == top_;
    }


    template<typename Value,
             typename Size>
    Size
    NODE_STACK::size() const
    {
      return top_ - bot_;
    }


    template<typename Value,
             typename Size>
    Size
    NODE_STACK::capacity() const
    {
      return stk_.size();
    }


  }

}


#undef NODE_STACK


#endif
