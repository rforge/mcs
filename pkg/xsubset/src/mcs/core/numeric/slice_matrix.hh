/**
 * @file slice_matrix.hh
 */

#ifndef MCS_CORE_NUMERIC_SLICE_MATRIX_HH
#define MCS_CORE_NUMERIC_SLICE_MATRIX_HH


#include "matrix_base.hh"


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
               typename Size>
      class matrix;



      template<typename Value,
               typename Size>
      class slice_matrix : public matrix_base<Value, Size, slice_matrix>
      {

        friend class matrix_base<Value, Size, matrix>;


      private:

        slice_matrix() = delete;

        slice_matrix(Value* base,
                     Size ldim,
                     Size nrow,
                     Size ncol);

        slice_matrix(const slice_matrix<Value, Size>& x);

      public:
        
        ~slice_matrix();

        template<template<typename V,
                          typename S>
                 class D>
        void
        copy(const matrix_base<Value, Size, D>& x);
        
        template<template<typename V,
                          typename S>
                 class D>
        void
        swap(matrix_base<Value, Size, D>&& x);

        void
        fill(Value x);

        slice_matrix<Value, Size>&
        operator =(const matrix<Value, Size>& x);

        slice_matrix<Value, Size>&
        operator =(const slice_matrix<Value, Size>& x);

        slice_matrix<Value, Size>&
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
  swap(mcs::core::numeric::slice_matrix<Value, Size>& a,
       mcs::core::numeric::slice_matrix<Value, Size>& b);


}


#include "slice_matrix.cc"
#endif
