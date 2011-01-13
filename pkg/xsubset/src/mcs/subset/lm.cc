/**
 * @file lm.cc
 */
#ifndef MCS_SUBSET_LM_CC
#define MCS_SUBSET_LM_CC


#include "../core/numeric/slice.hh"
#include "../core/numeric/matrix.hh"
#include "../core/numeric/lapack.hh"

#include "lm.hh"


#define LM lm<Value,                            \
              Size>


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    LM::lm(const Size nobs,
           const Size nvar,
           const Value* ay) :
      nobs_(nobs),
      nvar_(nvar),
      rz_(nobs, nvar + 1, ay),
      rss_(0)
    {
      mcs::core::numeric::lapack::geqrf(rz_);
      rss_ = rz_(nvar_, nvar_);
    }


    template<typename Value,
             typename Size>
    LM::lm(const mcs::core::numeric::matrix<Value, Size>& ay) :
      nobs_(ay.nrow()),
      nvar_(ay.ncol() - 1),
      rz_(ay),
      rss_(0)
    {
      mcs::core::numeric::lapack::geqrf(rz_);
      rss_ = rz_(nvar_, nvar_);
    }


    template<typename Value,
             typename Size>
    Size
    LM::nobs() const
    {
      return nobs_;
    }


    template<typename Value,
             typename Size>
    Size
    LM::nvar() const
    {
      return nvar_;
    }


    template<typename Value,
             typename Size>
    mcs::core::numeric::matrix<Value, Size>
    LM::rz() const
    {
      return rz_(mcs::core::numeric::slice<Size>(0, nvar_ + 1),
                 mcs::core::numeric::slice<Size>(0, nvar_ + 1));
    }


    template<typename Value,
             typename Size>
    Value
    LM::rss() const
    {
      return rss_;
    }


    template<typename Value,
             typename Size>
    Value
    LM::rms() const
    {
      return rss_ / (nobs_ - nvar_);
    }


  }

}


#undef LM


#endif
