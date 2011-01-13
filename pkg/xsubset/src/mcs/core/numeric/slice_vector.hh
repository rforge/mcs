/**
 * @file slice_vector.hh
 */
#ifndef MCS_CORE_NUMERIC_SLICE_VECTOR_HH
#define MCS_CORE_NUMERIC_SLICE_VECTOR_HH


#include "vector_base.hh"


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
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      class matrix_base;



      template<typename Value,
               typename Size>
      class matrix;



      template<typename Value,
               typename Size>
      class slice_matrix;



      template<typename Value,
               typename Size>
      class slice_vector : public vector_base<Value, Size, slice_vector>
      {

        friend class vector_base<Value, Size, vector>;
        friend class vector_base<Value, Size, numeric::slice_vector>;

        friend class vector<Value, Size>;

        friend class matrix_base<Value, Size, matrix>;
        friend class matrix_base<Value, Size, numeric::slice_matrix>;


      private:

        Size inc_;


      private:

        slice_vector() = delete;

        slice_vector(Value* base,
                     Size inc,
                     Size len);

        slice_vector(const slice_vector<Value, Size>& x);

      public:

        ~slice_vector();

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
        
        slice_vector<Value, Size>&
        operator =(const slice_vector<Value, Size>& x);

        template<template<typename V,
                          typename S>
                 class D>
        slice_vector<Value, Size>&
        operator =(const vector_base<Value, Size, D>& x);

        slice_vector<Value, Size>&
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
  swap(mcs::core::numeric::slice_vector<Value, Size>& x,
       mcs::core::numeric::slice_vector<Value, Size>& y);


}


#include "slice_vector.cc"
#endif
