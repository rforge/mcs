/**
 * @file preorder.hh
 */
#ifndef MCS_SUBSET_PREORDER_HH
#define MCS_SUBSET_PREORDER_HH


#include <vector>

#include "../core/numeric/matrix.hh"
#include "../core/numeric/givens.hh"

#include "node.hh"


namespace mcs
{

  namespace subset
  {

    template<typename Value,
             typename Size>
    class preorder
    {

    private:

      Size prad_;

      mutable std::vector<mcs::core::numeric::givens<Value> > gwork_;

      mutable std::vector<Value> vwork_;

      mutable std::vector<Size> swork_;

      mutable mcs::core::numeric::matrix<Value, Size> rwork_;


    public:

      preorder(Size nvar,
               Size mark,
               Size prad);

      void
      apply(node<Value, Size>& x) const;

      void
      apply(Size n,
            std::vector<Size>& s,
            Size k,
            mcs::core::numeric::matrix<Value, Size>& rz) const;

    };

  }

}


#include "preorder.cc"
#endif
