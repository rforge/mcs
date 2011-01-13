/**
 * @file vector_base.hh
 */
#ifndef MCS_CORE_NUMERIC_VECTOR_BASE_HH
#define MCS_CORE_NUMERIC_VECTOR_BASE_HH


#include "slice.hh"


namespace mcs
{

  namespace core
  {

    namespace numeric
    {



      template<typename Value,
               typename Size>
      class vector;



      template<typename Value,
               typename Size>
      class slice_vector;



      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      class vector_base
      {

        friend class vector_base<Value, Size, vector>;
        friend class vector_base<Value, Size, slice_vector>;

        friend class vector<Value, Size>;
        friend class slice_vector<Value, Size>;


      protected:

        Value* base_;

        Size len_;


      protected:

        vector_base();

        vector_base(Value* base,
                    Size len);

        vector_base(const vector_base<Value, Size, Derived>& x);

        ~vector_base();

      public:

        Value*
        base();

        const Value*
        base() const;

        Size
        inc() const;

        Size
        len() const;

        Value&
        at(Size pos);

        const Value&
        at(Size pos) const;
        
        slice_vector<Value, Size>
        at(const slice<Size>& s);

        const slice_vector<Value, Size>
        at(const slice<Size>& s) const;

        template<template<typename V,
                          typename S>
                 class D>
        void
        copy(const vector_base<Value, Size, D>& x);
        
        template<template<typename V,
                          typename S>
                 class D>
        void
        swap(vector_base<Value, Size, D>&& x);

        void
        fill(Value x);

        vector_base<Value, Size, Derived>&
        operator =(const vector_base<Value, Size, Derived>& x);

        template<template<typename V,
                          typename S>
                 class D>
        vector_base<Value, Size, Derived>&
        operator =(const vector_base<Value, Size, D>& x);

        template<template<typename V,
                          typename S>
                 class D>
        vector_base<Value, Size, Derived>&
        operator =(vector_base<Value, Size, D>&& x);

        vector_base<Value, Size, Derived>&
        operator =(Value x);
      
        Value&
        operator ()(Size pos);

        const Value&
        operator ()(Size pos) const;

        slice_vector<Value, Size>
        operator ()(const slice<Size>& pos);

        const slice_vector<Value, Size>
        operator ()(const slice<Size>& pos) const;

      };


    }

  }

}


namespace std
{


  template<typename Value,
           typename Size,
           template<typename V,
                    typename S>
           class Derived>
  void
  swap(mcs::core::numeric::vector_base<Value, Size, Derived>&& x,
       mcs::core::numeric::vector_base<Value, Size, Derived>&& y);


}


#include "vector_base.cc"
#endif
