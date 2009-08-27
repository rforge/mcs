#ifndef _MCS_XUTIL_VECTOR_CC_
#define _MCS_XUTIL_VECTOR_CC_


#include <algorithm>

#include "assert.hh"

#include "range.hh"
#include "vector.hh"
#include "subvector.hh"



namespace MCS
{

  namespace xutil
  {


    template<typename DataType>
    vector<DataType>::vector() :
      _stor(0),
      _base(0), _inc(-1), _len(-1)
    {
    }


    template<typename DataType>
    vector<DataType>::vector(int len) :
      _stor(new data_type[len]),
      _base(_stor), _inc(1), _len(len)
    {
      MCS_ASSERT(len >= 0);
    }


    template<typename DataType>
    vector<DataType>::vector(int len, data_type x) :
      _stor(new data_type[len]),
      _base(_stor), _inc(1), _len(len)
    {
      MCS_ASSERT(len >= 0);

      copy(x);
    }


    template<typename DataType>
    vector<DataType>::vector(const_vector_ref v) :
      _stor(0),
      _base(0), _inc(-1), _len(-1)
    {
      if (v._base)
        {
          _stor = new data_type[v._len];
          _base = _stor;
          _inc = 1;
          _len = v._len;

          copy(v);
        }
    }


    template<typename DataType>
    vector<DataType>::vector(data_ptr stor, int len) :
      _stor(stor),
      _base(stor), _inc(1), _len(len)
    {
      MCS_ASSERT(stor);
      MCS_ASSERT(len >= 0);
    }


    template<typename DataType>
    vector<DataType>::vector(data_ptr base, int inc, int len) :
      _stor(0),
      _base(base), _inc(inc), _len(len)
    {
      MCS_ASSERT(base);
      MCS_ASSERT(inc > 0);
      MCS_ASSERT(len >= 0);
    }


    template<typename DataType>
    vector<DataType>::~vector()
    {
      if (_stor)
        {
          delete [] _stor;
        }
    }


    template<typename DataType>
    typename vector<DataType>::vector_ref
    vector<DataType>::operator =(data_type x)
    {
      copy(x);

      return *this;
    }


    template<typename DataType>
    typename vector<DataType>::vector_ref
    vector<DataType>::operator =(vector_type v)
    {
      move(v);

      return *this;
    }


    template<typename DataType>
    typename vector<DataType>::data_ref
    vector<DataType>::operator ()(int i)
    {
      MCS_ASSERT(i >= 0 && i < _len);

      return *(_base + i * _inc);
    }


    template<typename DataType>
    typename vector<DataType>::const_data_ref
    vector<DataType>::operator ()(int i) const
    {
      MCS_ASSERT(i >= 0 && i < _len);

      return *(_base + i * _inc);
    }


    template<typename DataType>
    typename vector<DataType>::subvector_ref
    vector<DataType>::operator ()(range r)
    {
      MCS_ASSERT(r._start <= _len && r._end <= _len);

      if (r._end == -1)
        {
          r._end = _len;
        }

      return *(new subvector_type(_base + r._start * _inc,
                                  _inc,
                                  r._end - r._start));
    }


    template<typename DataType>
    typename vector<DataType>::const_subvector_ref
    vector<DataType>::operator ()(range r) const
    {
      MCS_ASSERT(r._start <= _len && r._end <= _len);

      if (r._end == -1)
        {
          r._end = _len;
        }

      return *(new subvector_type(_base + r._start * _inc,
                                  _inc,
                                  r._end - r._start));
    }


    template<typename DataType>
    int
    vector<DataType>::len() const
    {
      return _len;
    }


    template<typename DataType>
    int
    vector<DataType>::inc() const
    {
      return _inc;
    }


    template<typename DataType>
    void
    vector<DataType>::swap(vector_ref v)
    {
      move(v);
    }


    template<typename DataType>
    void
    vector<DataType>::swap(subvector_ref v)
    {
      MCS_ASSERT(_len == v._len);

      data_ptr ptr1 = _base;
      data_ptr ptr2 = v._base;
      for (int i = 0; i < _len; ++i)
        {
          std::swap(*ptr1, *ptr2);
          ptr1 += _inc;
          ptr2 += v._inc;
        }
    }


    template<typename DataType>
    void
    vector<DataType>::copy(data_type x)
    {
      data_ptr dst = _base;
      for (int i = 0; i < _len; ++i)
        {
          *dst = x;
          dst += _inc;
        }
    }


    template<typename DataType>
    void
    vector<DataType>::copy(const_vector_ref v)
    {
      data_ptr dst = _base;
      data_ptr src = v._base;
      for (int i = 0; i < _len; ++i)
        {
          *dst = *src;
          dst += _inc;
          src += v._inc;
        }
    }


    template<typename DataType>
    void
    vector<DataType>::move(vector_ref v)
    {
      std::swap(_stor, v._stor);
      std::swap(_base, v._base);
      std::swap(_inc, v._inc);
      std::swap(_len, v._len);
    }


  }

}


#endif
