/**
 * @file slice_matrix.cc
 */
#ifndef MCS_CORE_NUMERIC_SLICE_MATRIX_CC
#define MCS_CORE_NUMERIC_SLICE_MATRIX_CC


#include <algorithm>

#include "../../mcs.hh"

#include "matrix_base.hh"
#include "slice_matrix.hh"


#define SLICE_MATRIX slice_matrix<Value,        \
                                  Size>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
               typename Size>
      SLICE_MATRIX::slice_matrix(Value* const base,
                                 const Size ldim,
                                 const Size nrow,
                                 const Size ncol) :
        matrix_base<Value, Size, numeric::slice_matrix>(base, ldim, nrow, ncol)
      {
      }


      template<typename Value,
               typename Size>
      SLICE_MATRIX::slice_matrix(const slice_matrix<Value, Size>& x) :
        matrix_base<Value, Size, numeric::slice_matrix>(x)
      {
      }
      

      template<typename Value,
               typename Size>
      SLICE_MATRIX::~slice_matrix()
      {
        this->base_ = 0;
        this->ldim_ = 0;
        this->nrow_ = 0;
        this->ncol_ = 0;
      }


      template<typename Value,
               typename Size>
      template<template<typename V,
                        typename S>
               class D>
      void
      SLICE_MATRIX::copy(const matrix_base<Value, Size, D>& x)
      {
        MCS_ASSERT(x.nrow_ == this->nrow_, "invalid argument");
        MCS_ASSERT(x.ncol_ == this->ncol_, "invalid argument");

        const Value* src = x.base_;
        Value* dst = this->base_;
        const Value* last = this->base_ + this->ncol_ * this->ldim_;

        while (dst != last)
          {
            std::copy(src, src + this->nrow_, dst);
            src += x.ldim_;
            dst += this->ldim_;
          }
      }


      template<typename Value,
               typename Size>
      template<template<typename V,
                        typename S>
               class D>
      void
      SLICE_MATRIX::swap(matrix_base<Value, Size, D>&& x)
      {
        MCS_ASSERT(x.nrow_ == this->nrow_, "invalid argument");
        MCS_ASSERT(x.ncol_ == this->ncol_, "invalid argument");

        const Value* ptr1 = x.base_;
        Value* ptr2 = this->base_;
        Value* const last = this->base_ + this->ncol_ * this->ldim_;

        while (ptr2 != last)
          {
            std::swap_ranges(ptr1, ptr1 + this->nrow_, ptr2);
            ptr1 += x.ldim_;
            ptr2 += this->ldim_;
          }
      }


      template<typename Value,
               typename Size>
      void
      SLICE_MATRIX::fill(const Value x)
      {
        Value* dst = this->base_;
        Value* const last = this->base_ + this->ncol_ * this->ldim_;

        while (dst != last)
          {
            std::fill_n(dst, this->nrow_, x);
            dst += this->ldim_;
          }
      }


      template<typename Value,
               typename Size>
      slice_matrix<Value, Size>&
      SLICE_MATRIX::operator =(const matrix<Value, Size>& x)
      {
        copy(x);

        return *this;
      }


      template<typename Value,
               typename Size>
      slice_matrix<Value, Size>&
      SLICE_MATRIX::operator =(const slice_matrix<Value, Size>& x)
      {
        copy(x);

        return *this;
      }


      template<typename Value,
               typename Size>
      slice_matrix<Value, Size>&
      SLICE_MATRIX::operator =(const Value x)
      {
        fill(x);

        return *this;
      }


    }

  }

}


#undef SLICE_MATRIX


namespace std
{


  template<typename Value,
           typename Size>
  void
  swap(mcs::core::numeric::slice_matrix<Value, Size>& a,
       mcs::core::numeric::slice_matrix<Value, Size>& b)
  {
    a.swap(b);
  }


}


#endif
