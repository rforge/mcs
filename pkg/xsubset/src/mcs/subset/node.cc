/**
 * @file node.cc
 */
#ifndef MCS_SUBSET_NODE_CC
#define MCS_SUBSET_NODE_CC


#include <algorithm>

#include "../mcs.hh"

#include "../core/numeric/matrix.hh"
#include "../core/numeric/slice.hh"
#include "../core/numeric/givens.hh"
#include "../core/util/math.hh"
#include "../core/util/algo.hh"

#include "lm.hh"
#include "node.hh"


#define NODE node<Value, Size>


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    NODE::node(const Size nvar) :
      nvar_(nvar),
      s_(nvar, 0),
      mark_(0),
      rz_(nvar + 1, nvar + 1)
    {
    }


    template<typename Value,
             typename Size>
    NODE::node(const lm<Value, Size>& x,
               const Size mark) :
      nvar_(x.nvar()),
      s_(nvar_),
      mark_(mark),
      rz_(x.rz())
    {
      mcs::core::util::algo::iota_n(s_.begin(), nvar_, 0);
    }


    template<typename Value,
             typename Size>
    Size
    NODE::nvar() const
    {
      return nvar_;
    }


    template<typename Value,
             typename Size>
    Size
    NODE::mark() const
    {
      return mark_;
    }


    template<typename Value,
             typename Size>
    Value
    NODE::bound() const
    {
      return mcs::core::util::math::sqr(rz_(nvar_, nvar_));
    }


    template<typename Value,
             typename Size>
    void
    NODE::drop(const Size j,
               node<Value, Size>& x) const
    {
      MCS_ASSERT(j >= mark_, "invalid argument");
      MCS_ASSERT(j < nvar_, "invalid argument");

      drop(j, x.s_, x.rz_);
      x.nvar_ = nvar_ - 1;
      x.mark_ = j;
    }


    template<typename Value,
             typename Size>
    void
    NODE::drop(Size j,
               std::vector<Size>& s,
               mcs::core::numeric::matrix<Value, Size>& rz) const
    {
      // drop variable index
      std::copy_n(s_.begin(), j, s.begin());
      std::copy_n(s_.begin() + j + 1, nvar_ - j - 1, s.begin() + j);

      // drop column and restore QRD
      mcs::core::numeric::givens<Value>::zero(rz_(j    , mcs::core::numeric::slice<Size>(j + 1, nvar_ + 1)),
                                              rz_(j + 1, mcs::core::numeric::slice<Size>(j + 1, nvar_ + 1)),
                                              rz (j    , mcs::core::numeric::slice<Size>(j    , nvar_    )),
                                              rz (j + 1, mcs::core::numeric::slice<Size>(j    , nvar_    )));

      for (++j; j < nvar_; ++j)
	{
          mcs::core::numeric::givens<Value>::zero(rz (j    , mcs::core::numeric::slice<Size>(j    , nvar_    )),
                                                  rz_(j + 1, mcs::core::numeric::slice<Size>(j + 1, nvar_ + 1)),
                                                  rz (j    , mcs::core::numeric::slice<Size>(j    , nvar_    )),
                                                  rz (j + 1, mcs::core::numeric::slice<Size>(j    , nvar_    )));
	}
    }


    template<typename Value,
             typename Size>
    void
    NODE::swap(node<Value, Size>& x)
    {
      std::swap(nvar_, x.nvar_);
      std::swap(s_, x.s_);
      std::swap(mark_, x.mark_);
      std::swap(rz_, x.rz_);
    }


  }

}


#undef NODE


#endif
