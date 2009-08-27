#ifndef _MCS_XUTIL_SUBMATRIX_CC_
#define _MCS_XUTIL_SUBMATRIX_CC_


#include "matrix.hh"
#include "submatrix.hh"


#define MCS_XUTIL_SUBMATRIX_POOL_SIZE 1024


namespace MCS
{

  namespace xutil
  {


    template<typename DataType>
    struct __submatrix_pool
    {
      static const int SIZE = MCS_XUTIL_SUBMATRIX_POOL_SIZE;

      static submatrix<DataType> buffer[];

      static submatrix<DataType>* begin;
      static submatrix<DataType>* end;
      static submatrix<DataType>* current;

      static void* new_submatrix()
      {
        ++current;
        if (current == end)
          {
            current = begin;
          }
        return current;
      }
    };


    template<typename DataType>
    submatrix<DataType>
    __submatrix_pool<DataType>::buffer[__submatrix_pool<DataType>::SIZE];


    template<typename DataType>
    submatrix<DataType>*
    __submatrix_pool<DataType>::begin = __submatrix_pool<DataType>::buffer;


    template<typename DataType>
    submatrix<DataType>*
    __submatrix_pool<DataType>::end = __submatrix_pool<DataType>::begin + __submatrix_pool<DataType>::SIZE;


    template<typename DataType>
    submatrix<DataType>*
    __submatrix_pool<DataType>::current = __submatrix_pool<DataType>::end - 1;






    template<typename DataType>
    submatrix<DataType>::submatrix()
      : matrix<DataType>()
    {
    }


    template<typename DataType>
    submatrix<DataType>::submatrix(data_ptr base, int ldim, int nrow, int ncol)
      : matrix<DataType>(base, ldim, nrow, ncol)
    {
    }


    template<typename DataType>
    submatrix<DataType>::submatrix(const_submatrix_ref m)
      : matrix<DataType>(m)
    {
    }


    template<typename DataType>
    submatrix<DataType>::~submatrix()
    {
    }


    template<typename DataType>
    void*
    submatrix<DataType>::operator new(size_t)
    {
      return __submatrix_pool<DataType>::new_submatrix();
    }


    template<typename DataType>
    typename submatrix<DataType>::submatrix_ref
    submatrix<DataType>::operator =(data_type x)
    {
      matrix_type::copy(x);

      return *this;
    }


    template<typename DataType>
    typename submatrix<DataType>::submatrix_ref
    submatrix<DataType>::operator =(const_matrix_ref m)
    {
      matrix_type::copy(m);

      return *this;
    }


    template<typename DataType>
    typename submatrix<DataType>::submatrix_ref
    submatrix<DataType>::operator =(const_submatrix_ref m)
    {
      matrix_type::copy(m);

      return *this;
    }


    template<typename DataType>
    void
    submatrix<DataType>::swap(matrix_ref m)
    {
      MCS_ASSERT(matrix_type::_nrow == m._nrow);
      MCS_ASSERT(matrix_type::_ncol == m._ncol);

      data_ptr ptr1 = matrix_type::_base;
      data_ptr ptr2 = m._base;
      for (int j = 0; j < matrix_type::_ncol; ++j)
        {
          for (int i = 0; i < matrix_type::_nrow; ++i)
            {
              std::swap(ptr1[i], ptr2[i]);
            }
          ptr1 += matrix_type::_ldim;
          ptr2 += m._ldim;
        }
    }


  }

}


#endif
