/**
 * @file preorder.cc
 */
#ifndef MCS_SUBSET_PREORDER_CC
#define MCS_SUBSET_PREORDER_CC


#include <algorithm>
#include <functional>
#include <cstdlib>

#include "../core/numeric/matrix.hh"
#include "../core/numeric/blas.hh"
#include "../core/numeric/lapack.hh"
#include "../core/util/algo.hh"

#include "node.hh"
#include "preorder.hh"


#define PREORDER preorder<Value, Size>


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    PREORDER::preorder(const Size nvar,
                       const Size mark,
                       const Size prad) :
      prad_(0),
      gwork_(nvar + 1),
      vwork_(nvar + 1),
      swork_(nvar + 1),
      rwork_(nvar + 1, nvar + 1)
    {
      if (prad < (nvar - mark))
        prad_ = (nvar - mark) - prad;
    }


    template<typename Value,
             typename Size>
    void
    PREORDER::apply(node<Value, Size>& x) const
    {
      if ((x.nvar_ - x.mark_) <= prad_)
        return;
      apply(x.nvar_, x.s_, x.mark_, x.rz_);
    }


    template<typename Value,
             typename Size>
    void
    PREORDER::apply(const Size n,
                    std::vector<Size>& s,
                    const Size k,
                    mcs::core::numeric::matrix<Value, Size>& rz) const
    {
      // compute bounds
      for (Size j = k; j < n; ++j)
        {
          for (Size i = j + 1; i <= n; ++i)
            {
              Value t = rz(j, i);
              for (Size g = j + 1; g < i; ++g)
                t = -gwork_[g].s() * t + gwork_[g].c() * rz(g, i);
              gwork_[i].gen(t, rz(i, i));
            }
          vwork_[j - k] = std::abs(gwork_[n].r());
        }

      // compute order
      mcs::core::util::algo::iota_n(swork_.begin(), n - k, 0);
      mcs::core::util::algo::order_n(swork_.begin(), n - k, vwork_.begin(), std::greater<Value>());

      // permute index
      mcs::core::util::algo::permute_n(s.begin() + k, n - k, swork_.begin());

      // permute columns
      std::swap(rz, rwork_);
      for (Size i = k; i < n; ++i)
      	{
      	  Size j = k + swork_[i - k];
          mcs::core::numeric::blas::copy(mcs::core::numeric::slice<Size>(k, j + 1), rwork_.col(j), rz.col(i));
      	  rz(mcs::core::numeric::slice<Size>(j + 1, n + 1), i).fill(0);
      	}
      mcs::core::numeric::blas::copy(mcs::core::numeric::slice<Size>(k, n + 1), rwork_.col(n), rz.col(n));

      // restore QRD
      mcs::core::numeric::lapack::geqrf(mcs::core::numeric::slice<Size>(k, n + 1), rz);
    }


  }

}


#undef PREORDER


#endif
