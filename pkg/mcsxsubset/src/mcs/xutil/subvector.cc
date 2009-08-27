#ifndef _MCS_XUTIL_SUBVECTOR_CC_
#define _MCS_XUTIL_SUBVECTOR_CC_


#include "vector.hh"
#include "subvector.hh"


#define MCS_XUTIL_SUBVECTOR_POOL_SIZE 1024


namespace MCS
{

  namespace xutil
  {


    template<typename DataType>
    struct __subvector_pool
    {
      static const int SIZE = MCS_XUTIL_SUBVECTOR_POOL_SIZE;

      static subvector<DataType> buffer[];

      static subvector<DataType>* begin;
      static subvector<DataType>* end;
      static subvector<DataType>* current;

      static void* new_subvector()
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
    subvector<DataType>
    __subvector_pool<DataType>::buffer[__subvector_pool<DataType>::SIZE];


    template<typename DataType>
    subvector<DataType>*
    __subvector_pool<DataType>::begin = __subvector_pool<DataType>::buffer;


    template<typename DataType>
    subvector<DataType>*
    __subvector_pool<DataType>::end = __subvector_pool<DataType>::begin + __subvector_pool<DataType>::SIZE;


    template<typename DataType>
    subvector<DataType>*
    __subvector_pool<DataType>::current = __subvector_pool<DataType>::end - 1;




    template<typename DataType>
    subvector<DataType>::subvector()
      : vector<DataType>()
    {
    }


    template<typename DataType>
    subvector<DataType>::subvector(data_ptr base, int inc, int len)
      : vector<DataType>(base, inc, len)
    {
    }


    template<typename DataType>
    subvector<DataType>::subvector(const_subvector_ref s)
      : vector<DataType>(s)
    {
    }


    template<typename DataType>
    subvector<DataType>::~subvector()
    {
    }


    template<typename DataType>
    void*
    subvector<DataType>::operator new(size_t)
    {
      return __subvector_pool<DataType>::new_subvector();
    }


    template<typename DataType>
    typename subvector<DataType>::subvector_ref
    subvector<DataType>::operator =(data_type x)
    {
      vector_type::copy(x);

      return *this;
    }


    template<typename DataType>
    typename subvector<DataType>::subvector_ref
    subvector<DataType>::operator =(const_vector_ref v)
    {
      vector_type::copy(v);

      return *this;
    }


    template<typename DataType>
    typename subvector<DataType>::subvector_ref
    subvector<DataType>::operator =(const_subvector_ref s)
    {
      vector_type::copy(s);

      return *this;
    }


    template<typename DataType>
    void
    subvector<DataType>::swap(vector_ref v)
    {
      MCS_ASSERT(vector_type::_len == v._len);

      data_ptr ptr1 = vector_type::_base;
      data_ptr ptr2 = v._base;
      for (int i = 0; i < vector_type::_len; ++i)
        {
          std::swap(*ptr1, *ptr2);
          ptr1 += vector_type::_inc;
          ptr2 += v._inc;
        }
    }


  }

}


#endif
