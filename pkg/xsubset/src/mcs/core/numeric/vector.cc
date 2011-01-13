/**
 * @file vector.cc
 */
#ifndef MCS_CORE_NUMERIC_VECTOR_CC
#define MCS_CORE_NUMERIC_VECTOR_CC


#include <algorithm>
#include <utility>

#include "../../mcs.hh"

#include "slice.hh"
#include "vector_base.hh"
#include "vector.hh"
#include "slice_vector.hh"


#define VECTOR vector<Value,                    \
                      Size>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
               typename Size>
      VECTOR::vector(const Size len) :
        vector_base<Value, Size, numeric::vector>(new Value[len], len)
      {
        //std::cout << "vector::vector(Size)" << std::endl;

        MCS_ASSERT(len >= 0, "invalid argument");
      }


      template<typename Value,
               typename Size>
      VECTOR::vector(const Size len,
                     const Value x) :
        vector_base<Value, Size, numeric::vector>(new Value[len], len)
      {
        //std::cout << "vector::vector(Size, Value)" << std::endl;

        MCS_ASSERT(len >= 0, "invalid argument");

        fill(x);
      }


      template<typename Value,
               typename Size>
      VECTOR::vector(const vector<Value, Size>& x) :
        vector_base<Value, Size, numeric::vector>(new Value[x.len_], x.len_)
      {
        //std::cout << "vector::vector(const vector&)" << std::endl;

        copy(x);
      }


      template<typename Value,
               typename Size>
      VECTOR::vector(vector<Value, Size>&& x) :
        vector_base<Value, Size, numeric::vector>()
      {
        //std::cout << "vector::vector(vector&&)" << std::endl;

        swap(std::move(x));
      }


      template<typename Value,
               typename Size>
      VECTOR::vector(vector_base<Value, Size, numeric::vector>&& x) :
        vector_base<Value, Size, numeric::vector>()
      {
        //std::cout << "vector::vector(vector_base&&)" << std::endl;

        swap(std::move(x));
      }

                     
      template<typename Value,
               typename Size>
      template<template<typename V,
                        typename S>
               class D>
      VECTOR::vector(const vector_base<Value, Size, D>& x) :
        vector_base<Value, Size, numeric::vector>(new Value[x.len_], x.len_)
      {
        //std::cout << "vector::vector(const vector_base&)" << std::endl;

        copy(x);
      }
      

      template<typename Value,
               typename Size>
      VECTOR::~vector()
      {
        //std::cout << "vector::~vector()" << std::endl;

        if (this->base_)
          {
            delete [] this->base_;
            this->base_ = 0;
          }
      }


      template<typename Value,
               typename Size>
      Size
      VECTOR::inc() const
      {
        return 1;
      }


      template<typename Value,
               typename Size>
      Value&
      VECTOR::at(const Size pos)
      {
        MCS_ASSERT((pos >= 0) && (pos < this->len_), "index out of range");

        return this->base_[pos];
      }


      template<typename Value,
               typename Size>
      const Value&
      VECTOR::at(const Size pos) const
      {
        MCS_ASSERT((pos >= 0) && (pos < this->len_), "index out of range");

        return this->base_[pos];
      }


      template<typename Value,
               typename Size>
      slice_vector<Value, Size>
      VECTOR::at(const slice<Size>& s)
      {
        MCS_ASSERT(s.pos_ + s.len_ <= this->len_, "index out of range");

        return slice_vector<Value, Size>(this->base_ + s.pos_, 1, s.len_);
      }


      template<typename Value,
               typename Size>
      const slice_vector<Value, Size>
      VECTOR::at(const slice<Size>& s) const
      {
        MCS_ASSERT(s.pos_ + s.len_ <= this->len_, "index out of range");

        return slice_vector<Value, Size>(this->base_ + s.pos_, 1, s.len_);
      }


      template<typename Value,
               typename Size>
      void
      VECTOR::copy(const vector<Value, Size>& x)
      {
        //std::cout << "vector::copy(const vector&)" << std::endl;

        MCS_ASSERT(this->len_ == x.len_, "invalid argument");

        std::copy(x.base_, x.base_ + this->len_, this->base_);
      }


      template<typename Value,
               typename Size>
      void
      VECTOR::copy(const slice_vector<Value, Size>& x)
      {
        //std::cout << "vector::copy(const slice_vector&)" << std::endl;

        MCS_ASSERT(this->len_ == x.len_, "invalid argument");

        const Value* src = x.base_;
        Value* dst = this->base_;
        const Value* last = this->base_ + this->len_;

        while (dst != last)
          {
            *dst = *src;
            src += x.inc_;
            ++dst;
          }
      }


      template<typename Value,
               typename Size>
      template<template<typename V,
                        typename S>
               class D>
      void
      VECTOR::copy(const vector_base<Value, Size, D>& x)
      {
        //std::cout << "vector::copy(const vector_base&)" << std::endl;

        copy(*static_cast<const D<Value, Size>*>(&x));
      }


      template<typename Value,
               typename Size>
      void
      VECTOR::swap(vector<Value, Size>&& x)
      {
        //std::cout << "vector::swap(vector&&)" << std::endl;

        std::swap(this->base_, x.base_);
        std::swap(this->len_, x.len_);
      }


      template<typename Value,
               typename Size>
      void
      VECTOR::swap(slice_vector<Value, Size>&& x)
      {
        //std::cout << "vector::swap(slice_vector&&)" << std::endl;

        MCS_ASSERT(this->len_ == x.len_, "invalid argument");

        Value* ptr1 = x.base_;
        Value* ptr2 = this->base_;
        const Value* last = this->base_ + this->len_;

        while (ptr2 != last)
          {
            std::swap(*ptr1, *ptr2);
            ptr1 += x.inc_;
            ++ptr2;
          }
      }


      template<typename Value,
               typename Size>
      template<template<typename V,
                        typename S>
               class D>
      void
      VECTOR::swap(vector_base<Value, Size, D>&& x)
      {
        //std::cout << "vector::swap(vector_base&&)" << std::endl;

        swap(std::move(*static_cast<const D<Value, Size>*>(&x)));
      }


      template<typename Value,
               typename Size>
      void
      VECTOR::fill(const Value x)
      {
        //std::cout << "vector::fill(Value)" << std::endl;

        std::fill_n(this->base_, this->len_, x);
      }


      template<typename Value,
               typename Size>
      vector<Value, Size>&
      VECTOR::operator =(vector<Value, Size> x)
      {
        //std::cout << "vector::operator =(vector)" << std::endl;

        MCS_ASSERT(this->len_ == x.len_, "invalid argument");

        swap(x);

        return *this;
      }


      template<typename Value,
               typename Size>
      vector<Value, Size>&
      VECTOR::operator =(const Value x)
      {
        //std::cout << "vector::operator =(Value)" << std::endl;

        fill(x);

        return *this;
      }


    }

  }

}


#undef VECTOR


namespace std
{


  template<typename Value,
           typename Size>
  void
  swap(mcs::core::numeric::vector<Value, Size>& x,
       mcs::core::numeric::vector<Value, Size>& y)
  {
    x.swap(y);
  }


}


#endif
