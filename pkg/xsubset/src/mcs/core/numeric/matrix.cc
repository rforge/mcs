/**
 * @file matrix.cc
 */
#ifndef MCS_CORE_NUMERIC_MATRIX_CC
#define MCS_CORE_NUMERIC_MATRIX_CC


#include <algorithm>
#include <utility>

#include "../../mcs.hh"

#include "slice.hh"
#include "matrix_base.hh"
#include "matrix.hh"
#include "slice_matrix.hh"


#define MATRIX matrix<Value, Size>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
               typename Size>
      MATRIX::matrix(const Size nrow,
                     const Size ncol) :
        matrix_base<Value, Size, numeric::matrix>(new Value[nrow * ncol],
                                                  nrow, nrow, ncol)
      {
        MCS_ASSERT(nrow >= 0, "invalid argument");
        MCS_ASSERT(ncol >= 0, "invalid argument");
      }


      template<typename Value,
               typename Size>
      MATRIX::matrix(const Size nrow,
                     const Size ncol,
                     const Value x) :
        matrix_base<Value, Size, numeric::matrix>(new Value[nrow * ncol],
                                                  nrow, nrow, ncol)
      {
        MCS_ASSERT(nrow >= 0, "invalid argument");
        MCS_ASSERT(ncol >= 0, "invalid argument");

        fill(x);
      }


      template<typename Value,
               typename Size>
      MATRIX::matrix(const Size nrow,
                     const Size ncol,
                     const Value* const x) :
        matrix_base<Value, Size, numeric::matrix>(new Value[nrow * ncol],
                                                  nrow, nrow, ncol)
      {
        MCS_ASSERT(nrow >= 0, "invalid argument");
        MCS_ASSERT(ncol >= 0, "invalid argument");

        std::copy(x, x + nrow * ncol, this->base_);
      }


      template<typename Value,
               typename Size>
      MATRIX::matrix(const matrix<Value, Size>& x) :
        matrix_base<Value, Size, numeric::matrix>(new Value[x.nrow_ * x.ncol_],
                                                  x.nrow_, x.nrow_, x.ncol_)
      {
        copy(x);
      }


      template<typename Value,
               typename Size>
      MATRIX::matrix(matrix<Value, Size>&& x) :
        matrix_base<Value, Size, numeric::matrix>()
      {
        swap(std::move(x));
      }

      
      template<typename Value,
               typename Size>
      MATRIX::matrix(matrix_base<Value, Size, numeric::matrix>&& x) :
        matrix_base<Value, Size, numeric::matrix>()
      {
        swap(std::move(x));
      }
      
      
      template<typename Value,
               typename Size>
      template<template<typename V,
                        typename S>
               class D>
      MATRIX::matrix(const matrix_base<Value, Size, D>& x) :
        matrix_base<Value, Size, numeric::matrix>(new Value[x.nrow_ * x.ncol_],
                                                  x.nrow_, x.nrow_, x.ncol_)
      {
        copy(x);
      }
      

      template<typename Value,
               typename Size>
      MATRIX::~matrix()
      {
        if (this->base_)
          {
            delete [] this->base_;
            this->base_ = 0;
          }
      }


      template<typename Value,
               typename Size>
      void
      MATRIX::copy(const matrix<Value, Size>& x)
      {
        MCS_ASSERT(x.nrow_ == this->nrow_, "invalid argument");
        MCS_ASSERT(x.ncol_ == this->ncol_, "invalid argument");

        std::copy(x.base_, x.base_ + this->nrow_ * this->ncol_, this->base_);
      }


      template<typename Value,
               typename Size>
      void
      MATRIX::copy(const slice_matrix<Value, Size>& x)
      {
        MCS_ASSERT(x.nrow_ == this->nrow_, "invalid argument");
        MCS_ASSERT(x.ncol_ == this->ncol_, "invalid argument");

        const Value* src = x.base_;
        Value* dst = this->base_;
        const Value* last = this->base_ + this->nrow_ * this->ncol_;

        while (dst != last)
          {
            dst = std::copy(src, src + this->nrow_, dst);
            src += x.ldim_;
          }
      }


      template<typename Value,
               typename Size>
      template<template<typename V,
                        typename S>
               class D>
      void
      MATRIX::copy(const matrix_base<Value, Size, D>& x)
      {
        copy(*static_cast<const D<Value, Size>*>(&x));
      }


      template<typename Value,
               typename Size>
      void
      MATRIX::swap(matrix<Value, Size>&& x)
      {
	std::swap(this->base_, x.base_);
	std::swap(this->ldim_, x.ldim_);
	std::swap(this->nrow_, x.nrow_);
	std::swap(this->ncol_, x.ncol_);
      }


      template<typename Value,
               typename Size>
      void
      MATRIX::swap(slice_matrix<Value, Size>&& x)
      {
        MCS_ASSERT(this->nrow_ == x.nrow_, "invalid argument");
        MCS_ASSERT(this->ncol_ == x.ncol_, "invalid argument");

        const Value* ptr1 = x.base_;
        Value* ptr2 = this->base_;
        const Value* last = this->base_ + this->nrow_ * this->ncol_;

        while (ptr2 != last)
          {
            ptr2 = std::swap_ranges(ptr1, ptr1 + this->nrow_, ptr2);
            ptr1 += x.ldim_;
          }
      }


      template<typename Value,
               typename Size>
      template<template<typename V,
                        typename S>
               class D>
      void
      MATRIX::swap(matrix_base<Value, Size, D>&& x)
      {
        swap(std::move(*static_cast<const D<Value, Size>*>(&x)));
      }


      template<typename Value,
               typename Size>
      void
      MATRIX::fill(const Value x)
      {
        std::fill_n(this->base_, this->nrow_ * this->ncol_, x);
      }


      template<typename Value,
               typename Size>
      matrix<Value, Size>&
      MATRIX::operator =(matrix<Value, Size> x)
      {
        MCS_ASSERT(this->nrow_ == x.nrow_, "invalid argument");
        MCS_ASSERT(this->ncol_ == x.ncol_, "invalid argument");

        swap(x);

        return *this;
      }


      template<typename Value,
               typename Size>
      matrix<Value, Size>&
      MATRIX::operator =(const Value x)
      {
        fill(x);

        return *this;
      }


    }

  }

}


#undef MATRIX


namespace std
{


  template<typename Value,
           typename Size>
  void
  swap(mcs::core::numeric::matrix<Value, Size>& a,
       mcs::core::numeric::matrix<Value, Size>& b)
  {
    a.swap(b);
  }


}


#endif
