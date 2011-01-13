/**
 * @file matrix.hh
 */
#ifndef MCS_CORE_NUMERIC_MATRIX_HH
#define MCS_CORE_NUMERIC_MATRIX_HH


#include "slice.hh"
#include "matrix_base.hh"


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
               typename Size>
      class matrix : public matrix_base<Value, Size, matrix>
      {

      public:

        matrix() = delete;

        matrix(Size nrow,
               Size ncol);

        matrix(Size nrow,
               Size ncol,
               Value x);

        matrix(Size nrow,
               Size ncol,
               const Value* x);

        matrix(const matrix<Value, Size>& x);

        matrix(matrix<Value, Size>&& x);

        matrix(matrix_base<Value, Size, numeric::matrix>&& x);

        template<template<typename V,
                          typename S>
                 class D>
        matrix(const matrix_base<Value, Size, D>& x);

        ~matrix();

        void
        copy(const matrix<Value, Size>& x);

        void
        copy(const slice_matrix<Value, Size>& x);

        template<template<typename V,
                          typename S>
                 class D>
        void
        copy(const matrix_base<Value, Size, D>& x);

        void
        swap(matrix<Value, Size>&& x);

        void
        swap(slice_matrix<Value, Size>&& x);

        template<template<typename V,
                          typename S>
                 class D>
        void
        swap(matrix_base<Value, Size, D>&& x);

        void
        fill(Value x);

        matrix<Value, Size>&
        operator =(matrix<Value, Size> x);

        matrix<Value, Size>&
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
  swap(mcs::core::numeric::matrix<Value, Size>& a,
       mcs::core::numeric::matrix<Value, Size>& b);


}


#include "matrix.cc"
#endif
