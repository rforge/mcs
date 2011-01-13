/**
 * @file vector.hh
 */
#ifndef MCS_CORE_NUMERIC_VECTOR_HH
#define MCS_CORE_NUMERIC_VECTOR_HH


#include "slice.hh"
#include "vector_base.hh"


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
               typename Size>
      class vector : public vector_base<Value, Size, numeric::vector>
      {

      public:

        vector() = delete;

        explicit
        vector(Size len);

        vector(Size len,
               Value x);

        vector(const vector<Value, Size>& x);

        vector(vector<Value, Size>&& x);

        vector(vector_base<Value, Size, numeric::vector>&& x);

        template<template<typename V,
                          typename S>
                 class D>
        vector(const vector_base<Value, Size, D>& x);
        
        ~vector();

        Size
        inc() const;

        Value&
        at(Size pos);

        const Value&
        at(Size pos) const;
        
        slice_vector<Value, Size>
        at(const slice<Size>& s);

        const slice_vector<Value, Size>
        at(const slice<Size>& s) const;

        void
        copy(const vector<Value, Size>& x);

        void
        copy(const slice_vector<Value, Size>& x);

        template<template<typename V,
                          typename S>
                 class D>
        void
        copy(const vector_base<Value, Size, D>& x);

        void
        swap(vector<Value, Size>&& x);

        void
        swap(slice_vector<Value, Size>&& x);

        template<template<typename V,
                          typename S>
                 class D>
        void
        swap(vector_base<Value, Size, D>&& x);

        void
        fill(Value x);

        vector<Value, Size>&
        operator =(vector<Value, Size> x);

        vector<Value, Size>&
        operator =(Value x);

      };


    }

  }

}


namespace std
{


  template<typename Value,
           typename Size>
  void
  swap(mcs::core::numeric::vector<Value, Size>& x,
       mcs::core::numeric::vector<Value, Size>& y);


}


#include "vector.cc"
#endif
