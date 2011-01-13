/**
 * @file slice_vector.cc
 */
#ifndef MCS_CORE_NUMERIC_SLICE_VECTOR_CC
#define MCS_CORE_NUMERIC_SLICE_VECTOR_CC


#include <utility>

#include "../../mcs.hh"

#include "vector_base.hh"
#include "slice_vector.hh"


#define SLICE_VECTOR slice_vector<Value,        \
                                  Size>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
               typename Size>
      SLICE_VECTOR::slice_vector(Value* const base,
                                 const Size inc,
                                 const Size len) :
        vector_base<Value, Size, numeric::slice_vector>(base, len),
        inc_(inc)
      {
        //std::cout << "slice_vector::slice_vector(Value*,Size,Size)" << std::endl;
      }


      template<typename Value,
               typename Size>
      SLICE_VECTOR::slice_vector(const slice_vector<Value, Size>& x) :
        vector_base<Value, Size, numeric::slice_vector>(x.base_, x.len_),
        inc_(x.inc_)
      {
        //std::cout << "slice_vector::slice_vector(const slice_vector&)" << std::endl;
      }
      

      template<typename Value,
               typename Size>
      SLICE_VECTOR::~slice_vector()
      {
        //std::cout << "slice_vector::~slice_vector()" << std::endl;

        this->base_ = 0;
        this->inc_ = 0;
        this->len_ = 0;
      }


      template<typename Value,
               typename Size>
      Size
      SLICE_VECTOR::inc() const
      {
        return inc_;
      }


      template<typename Value,
               typename Size>
      Value&
      SLICE_VECTOR::at(const Size pos)
      {
        MCS_ASSERT((pos >= 0) && (pos < this->len_), "index out of range");

        return this->base_[pos * this->inc_];
      }


      template<typename Value,
               typename Size>
      const Value&
      SLICE_VECTOR::at(const Size pos) const
      {
        MCS_ASSERT((pos >= 0) && (pos < this->len_), "index out of range");

        return this->base_[pos * inc_];
      }


      template<typename Value,
               typename Size>
      slice_vector<Value, Size>
      SLICE_VECTOR::at(const slice<Size>& s)
      {
        MCS_ASSERT(s.pos_ + s.len_ <= this->len_, "index out of range");

        return slice_vector<Value, Size>(this->base_ + s.pos_ * this->inc_,
                                         this->inc_, s.len_);
      }


      template<typename Value,
               typename Size>
      const slice_vector<Value, Size>
      SLICE_VECTOR::at(const slice<Size>& s) const
      {
        MCS_ASSERT(s.pos_ + s.len_ <= this->len_, "index out of range");

        return slice_vector<Value, Size>(this->base_ + s.pos_ * this->inc_,
                                         this->inc_, s.len_);
      }


      template<typename Value,
               typename Size>
      void
      SLICE_VECTOR::copy(const vector<Value, Size>& x)
      {
        //std::cout << "slice_vector::copy(const vector&)" << std::endl;

        MCS_ASSERT(this->len_ = x.len_, "invalid argument");

        Value* src = x.base_;
        const Value* last = x.base_ + x.len_;
        Value* dst = this->base_;

        while (src != last)
          {
            *dst = *src;
            ++src;
            dst += this->inc_;
          }
      }


      template<typename Value,
               typename Size>
      void
      SLICE_VECTOR::copy(const slice_vector<Value, Size>& x)
      {
        //std::cout << "slice_vector::copy(const slice_vector&)" << std::endl;

        MCS_ASSERT(this->len_ = x.len_, "invalid argument");

        Value* src = x.base_;
        Value* dst = this->base_;

        for (Size pos = 0; pos < this->len_; ++pos)
          {
            *dst = *src;
            src += x.inc_;
            dst += this->inc_;
          }
      }


      template<typename Value,
               typename Size>
      template<template<typename V,
                        typename S>
               class D>
      void
      SLICE_VECTOR::copy(const vector_base<Value, Size, D>& x)
      {
        //std::cout << "slice_vector::copy(const vector_base&)" << std::endl;

        copy(*static_cast<const D<Value, Size>*>(&x));
      }


      template<typename Value,
               typename Size>
      void
      SLICE_VECTOR::swap(vector<Value, Size>&& x)
      {
        //std::cout << "slice_vector::swap(vector&&)" << std::endl;

        MCS_ASSERT(this->len_ = x.len_, "invalid argument");

        Value* ptr1 = x.base_;
        const Value* last = x.base_ + x.len_;
        Value* ptr2 = this->base_;

        while (ptr1 != last)
          {
            std::swap(*ptr1, *ptr2);
            ++ptr1;
            ptr2 += this->inc_;
          }
      }


      template<typename Value,
               typename Size>
      void
      SLICE_VECTOR::swap(slice_vector<Value, Size>&& x)
      {
        //std::cout << "slice_vector::swap(slice_vector&&)" << std::endl;

        MCS_ASSERT(this->len_ = x.len_, "invalid argument");

        Value* ptr1 = x.base_;
        Value* ptr2 = this->base_;

        for (Size pos = 0; pos < this->len_; ++pos)
          {
            std::swap(*ptr1, *ptr2);
            ptr1 += x.inc_;
            ptr2 += this->inc_;
          }
      }


      template<typename Value,
               typename Size>
      template<template<typename V,
                        typename S>
               class D>
      void
      SLICE_VECTOR::swap(vector_base<Value, Size, D>&& x)
      {
        //std::cout << "slice_vector::swap(vector_base&&)" << std::endl;

        swap(std::move(*static_cast<const D<Value, Size>*>(&x)));
      }


      template<typename Value,
               typename Size>
      void
      SLICE_VECTOR::fill(const Value x)
      {
        //std::cout << "slice_vector::fill(Value)" << std::endl;

        Value* dst = this->base_;

        for (Size pos = 0; pos < this->len_; ++pos)
          {
            *(dst++) = x;
          }
      }


      template<typename Value,
               typename Size>
      slice_vector<Value, Size>&
      SLICE_VECTOR::operator =(const slice_vector<Value, Size>& x)
      {
        //std::cout << "slice_vector::operator=(const slice_vector&)" << std::endl;

        copy(x);

        return *this;
      }


      template<typename Value,
               typename Size>
      template<template<typename V,
                        typename S>
               class D>
      slice_vector<Value, Size>&
      SLICE_VECTOR::operator =(const vector_base<Value, Size, D>& x)
      {
        //std::cout << "slice_vector::operator=(const vector_base&)" << std::endl;

        copy(x);

        return *this;
      }


      template<typename Value,
               typename Size>
      slice_vector<Value, Size>&
      SLICE_VECTOR::operator =(const Value x)
      {
        //std::cout << "slice_vector::operator=(const Value)" << std::endl;

        fill(x);

        return *this;
      }


    }

  }

}


#undef SLICE_VECTOR


namespace std
{


  template<typename Value,
           typename Size>
  void
  swap(mcs::core::numeric::slice_vector<Value, Size>& x,
       mcs::core::numeric::slice_vector<Value, Size>& y)
  {
    x.swap(y);
  }


}


#endif
