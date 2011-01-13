/**
 * @file node.hh
 */
#ifndef MCS_SUBSET_NODE_HH
#define MCS_SUBSET_NODE_HH


#include <vector>

#include "../core/numeric/matrix.hh"

#include "lm.hh"
#include "criteria.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    class node_stack;



    template<typename Value,
             typename Size>
    class preorder;



    template<typename Value,
             typename Size>
    class subset_table;



    template<typename Value,
             typename Size,
             template<typename V,
                      typename S>
             class Criterion>
    class subset_table1;



    template<typename Value,
             typename Size>
    class node
    {

      friend class node_stack<Value, Size>;
      friend class preorder<Value, Size>;
      friend class subset_table<Value, Size>;

      friend class subset_table1<Value, Size, aic>;
      friend class subset_table1<Value, Size, cp>;


    private:

      Size nvar_;

      std::vector<Size> s_;

      Size mark_;

      mcs::core::numeric::matrix<Value, Size> rz_;


    public:

      node(Size nvar);

      node(const lm<Value, Size>& x,
           Size mark);


      Size
      nvar() const;

      Size
      mark() const;

      Value
      bound() const;

      void
      drop(Size j,
           node<Value, Size>& x) const;

      void
      drop(Size j,
           std::vector<Size>& s,
           mcs::core::numeric::matrix<Value, Size>& rz) const;

      void
      swap(node<Value, Size>& x);

    };


  }

}


#include "node.cc"
#endif
