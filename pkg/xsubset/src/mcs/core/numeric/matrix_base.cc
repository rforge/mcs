/**
 * @file matrix_base.cc
 */
#ifndef MCS_CORE_NUMERIC_MATRIX_BASE_CC
#define MCS_CORE_NUMERIC_MATRIX_BASE_CC


#include <utility>
#include <algorithm>

#include "../../mcs.hh"

#include "slice.hh"
#include "slice_vector.hh"
#include "matrix_base.hh"
#include "slice_matrix.hh"


#define MATRIX_BASE matrix_base<Value,                    \
                                Size,                     \
                                Derived>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      MATRIX_BASE::matrix_base():
	base_(0),
        ldim_(0),
	nrow_(0),
	ncol_(0)
      {
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      MATRIX_BASE::matrix_base(Value* const base,
                               const Size ldim,
                               const Size nrow,
                               const Size ncol) :
	base_(base),
        ldim_(ldim),
	nrow_(nrow),
	ncol_(ncol)
      {
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      MATRIX_BASE::matrix_base(const matrix_base<Value, Size, Derived>& x) :
	base_(x.base_),
        ldim_(x.ldim_),
	nrow_(x.nrow_),
	ncol_(x.ncol_)
      {
      }

 
      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      MATRIX_BASE::~matrix()
      {
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      Value*
      MATRIX_BASE::base()
      {
	return base_;
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const Value*
      MATRIX_BASE::base() const
      {
	return base_;
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      Size
      MATRIX_BASE::ldim() const
      {
        return ldim_;
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      Size
      MATRIX_BASE::nrow() const
      {
	return nrow_;
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      Size
      MATRIX_BASE::ncol() const
      {
	return ncol_;
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      Value&
      MATRIX_BASE::at(const Size irow,
                      const Size icol)
      {
        MCS_ASSERT((irow >= 0) && (irow < nrow_), "index out of range");
        MCS_ASSERT((icol >= 0) && (icol < ncol_), "index out of range");

	return base_[irow + icol * ldim_];
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const Value&
      MATRIX_BASE::at(const Size irow,
                      const Size icol) const
      {
        MCS_ASSERT((irow >= 0) && (irow < nrow_), "index out of range");
        MCS_ASSERT((icol >= 0) && (icol < ncol_), "index out of range");

	return base_[irow + icol * ldim_];
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_vector<Value, Size>
      MATRIX_BASE::at(const Size irow,
                      const numeric::slice<Size>& scol)
      {
	MCS_ASSERT((irow >= 0) && (irow < nrow_), "index out of range");
	MCS_ASSERT((scol.pos_ + scol.len_) <= ncol_, "index out of range");

	return slice_vector<Value, Size>(base_ + irow + scol.pos_ * ldim_,
                                         ldim_, scol.len_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_vector<Value, Size>
      MATRIX_BASE::at(const Size irow,
                      const numeric::slice<Size>& scol) const
      {
	MCS_ASSERT((irow >= 0) && (irow < nrow_), "index out of range");
	MCS_ASSERT((scol.pos_ + scol.len_) <= ncol_, "index out of range");

	return slice_vector<Value, Size>(base_ + irow + scol.pos_ * ldim_,
                                         ldim_, scol.len_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_vector<Value, Size>
      MATRIX_BASE::at(const numeric::slice<Size>& srow,
                      const Size icol)
      {
	MCS_ASSERT((srow.pos_ + srow.len_) <= nrow_, "index out of range");
	MCS_ASSERT((icol >= 0) && (icol < ncol_), "index out of range");

	return slice_vector<Value, Size>(base_ + srow.pos_ + icol * ldim_, 1,
                                         srow.len_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_vector<Value, Size>
      MATRIX_BASE::at(const numeric::slice<Size>& srow,
                      const Size icol) const
      {
	MCS_ASSERT((srow.pos_ + srow.len_) <= nrow_, "index out of range");
	MCS_ASSERT((icol >= 0) && (icol < ncol_), "index out of range");

	return slice_vector<Value, Size>(base_ + srow.pos_ + icol * ldim_, 1,
                                         srow.len_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_matrix<Value, Size>
      MATRIX_BASE::at(const numeric::slice<Size>& srow,
                      const numeric::slice<Size>& scol)
      {
	MCS_ASSERT((srow.pos_ + srow.len_) <= nrow_, "index out of range");
	MCS_ASSERT((scol.pos_ + scol.len_) <= ncol_, "index out of range");

	return slice_matrix<Value, Size>(base_ + srow.pos_ + scol.pos_ * ldim_,
                                         ldim_, srow.len_, scol.len_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_matrix<Value, Size>
      MATRIX_BASE::at(const numeric::slice<Size>& srow,
                      const numeric::slice<Size>& scol) const
      {
	MCS_ASSERT((srow.pos_ + srow.len_) <= nrow_, "index out of range");
	MCS_ASSERT((scol.pos_ + scol.len_) <= ncol_, "index out of range");

	return slice_matrix<Value, Size>(base_ + srow.pos_ + scol.pos_ * ldim_,
                                         ldim_, srow.len_, scol.len_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_vector<Value, Size>
      MATRIX_BASE::row(const Size irow)
      {
	MCS_ASSERT((irow >= 0) && (irow < nrow_), "index out of range");

	return slice_vector<Value, Size>(base_ + irow, ldim_, ncol_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_vector<Value, Size>
      MATRIX_BASE::row(const Size irow) const
      {
	MCS_ASSERT((irow >= 0) && (irow < nrow_), "index out of range");

	return slice_vector<Value, Size>(base_ + irow, ldim_, ncol_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_vector<Value, Size>
      MATRIX_BASE::col(const Size icol)
      {
	MCS_ASSERT((icol >= 0) && (icol < ncol_), "index out of range");

	return slice_vector<Value, Size>(base_ + icol * ldim_, 1, nrow_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_vector<Value, Size>
      MATRIX_BASE::col(const Size icol) const
      {
	MCS_ASSERT((icol >= 0) && (icol < ncol_), "index out of range");

	return slice_vector<Value, Size>(base_ + icol * ldim_, 1, nrow_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_matrix<Value, Size>
      MATRIX_BASE::rows(const numeric::slice<Size>& srow)
      {
	MCS_ASSERT((srow.pos_ + srow.len_) <= nrow_, "index out of range");

	return slice_matrix<Value, Size>(base_ + srow.pos_, ldim_, srow.len_,
                                         ncol_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_matrix<Value, Size>
      MATRIX_BASE::rows(const numeric::slice<Size>& srow) const
      {
	MCS_ASSERT((srow.pos_ + srow.len_) <= nrow_, "index out of range");

	return slice_matrix<Value, Size>(base_ + srow.pos_, ldim_, srow.len_,
                                         ncol_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_matrix<Value, Size>
      MATRIX_BASE::cols(const numeric::slice<Size>& scol)
      {
	MCS_ASSERT((scol.pos_ + scol.len_) <= ncol_, "index out of range");

	return slice_matrix<Value, Size>(base_ + scol.pos_ * ldim_, ldim_,
                                         nrow_, scol.len_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_matrix<Value, Size>
      MATRIX_BASE::cols(const numeric::slice<Size>& scol) const
      {
	MCS_ASSERT((scol.pos_ + scol.len_) <= ncol_, "index out of range");

	return slice_matrix<Value, Size>(base_ + scol.pos_ * ldim_, ldim_,
                                         nrow_, scol.len_);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_vector<Value, Size>
      MATRIX_BASE::diag()
      {
	return slice_vector<Value, Size>(base_, ldim_ + 1,
                                         std::min(nrow_, ncol_));
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_vector<Value, Size>
      MATRIX_BASE::diag() const
      {
	return slice_vector<Value, Size>(base_, ldim_ + 1,
                                         std::min(nrow_, ncol_));
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      void
      MATRIX_BASE::fill(const Value x)
      {
        return static_cast<Derived<Value, Size>*>(this)->fill(x);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      template<template<typename V,
                        typename S>
               class D>
      void
      MATRIX_BASE::copy(const matrix_base<Value, Size, D>& x)
      {
        return static_cast<Derived<Value, Size>*>(this)->copy(*static_cast<const D<Value, Size>*>(&x));
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      template<template<typename V,
                        typename S>
               class D>
      void
      MATRIX_BASE::swap(matrix_base<Value, Size, D>&& x)
      {
        return static_cast<Derived<Value, Size>*>(this)->swap(*static_cast<const D<Value, Size>*>(&x));
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      matrix_base<Value, Size, Derived>&
      MATRIX_BASE::operator =(const matrix_base<Value, Size, Derived>& x)
      {
        return static_cast<Derived<Value, Size>*>(this)->operator =(*static_cast<const Derived<Value, Size>*>(&x));
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      template<template<typename V,
                        typename S>
               class D>
      matrix_base<Value, Size, Derived>&
      MATRIX_BASE::operator =(const matrix_base<Value, Size, D>& x)
      {
        return static_cast<Derived<Value, Size>*>(this)->operator =(*static_cast<const D<Value, Size>*>(&x));
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      template<template<typename V,
                        typename S>
               class D>
      matrix_base<Value, Size, Derived>&
      MATRIX_BASE::operator =(matrix_base<Value, Size, D>&& x)
      {
        return static_cast<Derived<Value, Size>*>(this)->operator =(std::move(*static_cast<const D<Value, Size>*>(&x)));
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      matrix_base<Value, Size, Derived>&
      MATRIX_BASE::operator =(const Value x)
      {
        fill(x);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      Value&
      MATRIX_BASE::operator ()(const Size irow,
                               const Size icol)
      {
	return at(irow, icol);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const Value&
      MATRIX_BASE::operator ()(const Size irow,
                               const Size icol) const
      {
	return at(irow, icol);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_vector<Value, Size>
      MATRIX_BASE::operator ()(const Size irow,
                               const numeric::slice<Size>& scol)
      {
        return at(irow, scol);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_vector<Value, Size>
      MATRIX_BASE::operator ()(const Size irow,
                               const numeric::slice<Size>& scol) const
      {
        return at(irow, scol);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_vector<Value, Size>
      MATRIX_BASE::operator ()(const numeric::slice<Size>& srow,
                               const Size icol)
      {
        return at(srow, icol);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_vector<Value, Size>
      MATRIX_BASE::operator ()(const numeric::slice<Size>& srow,
                               const Size icol) const
      {
        return at(srow, icol);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_matrix<Value, Size>
      MATRIX_BASE::operator ()(const numeric::slice<Size>& srow,
                               const numeric::slice<Size>& scol)
      {
        return at(srow, scol);
      }


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_matrix<Value, Size>
      MATRIX_BASE::operator ()(const numeric::slice<Size>& srow,
                               const numeric::slice<Size>& scol) const
      {
        return at(srow, scol);
      }


    }

  }

}


#undef MATRIX_BASE


namespace std
{


  template<typename Value,
           typename Size,
           template<typename V,
                    typename S>
           class Derived>
  void
  swap(mcs::core::numeric::matrix_base<Value, Size, Derived>& a,
       mcs::core::numeric::matrix_base<Value, Size, Derived>& b)
  {
    a.swap(b);
  }


}


#endif
