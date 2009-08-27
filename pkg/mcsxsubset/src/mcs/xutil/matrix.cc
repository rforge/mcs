#ifndef _MCS_XUTIL_MATRIX_CC_
#define _MCS_XUTIL_MATRIX_CC_


#include <algorithm>

#include "assert.hh"

#include "range.hh"

#include "vector.hh"
#include "subvector.hh"

#include "matrix.hh"
#include "submatrix.hh"


namespace MCS
{

  namespace xutil
  {


    template<typename DataType>
    matrix<DataType>::matrix() :
      _stor(0),
      _base(0), _ldim(-1), _nrow(-1), _ncol(-1)
    {
    }


    template<typename DataType>
    matrix<DataType>::matrix(int nrow, int ncol) :
      _stor(new data_type[nrow * ncol]),
      _base(_stor), _ldim(nrow), _nrow(nrow), _ncol(ncol)
    {
      MCS_ASSERT(nrow >= 0);
      MCS_ASSERT(ncol >= 0);
    }


    template<typename DataType>
    matrix<DataType>::matrix(int nrow, int ncol, data_type x) :
      _stor(new data_type[nrow * ncol]),
      _base(_stor), _ldim(nrow), _nrow(nrow), _ncol(ncol)
    {
      MCS_ASSERT(nrow >= 0);
      MCS_ASSERT(ncol >= 0);

      copy(x);
    }


    template<typename DataType>
    matrix<DataType>::matrix(int nrow, int ncol, data_ptr x) :
      _stor(0),
      _base(x), _ldim(nrow), _nrow(nrow), _ncol(ncol)
    {
      MCS_ASSERT(x);
      MCS_ASSERT(nrow >= 0);
      MCS_ASSERT(ncol >= 0);
    }


    template<typename DataType>
    matrix<DataType>::matrix(data_ptr stor, int nrow, int ncol) :
      _stor(stor),
      _base(stor), _ldim(nrow), _nrow(nrow), _ncol(ncol)
    {
      MCS_ASSERT(stor);
      MCS_ASSERT(nrow >= 0);
      MCS_ASSERT(ncol >= 0);
    }


    template<typename DataType>
    matrix<DataType>::matrix(data_ptr base, int ldim, int nrow, int ncol) :
      _stor(0),
      _base(base), _ldim(ldim), _nrow(nrow), _ncol(ncol)
    {
      MCS_ASSERT(base);
      MCS_ASSERT(ldim > 0);
      MCS_ASSERT((nrow >= 0) && (nrow <= ldim));
      MCS_ASSERT(ncol >= 0);
    }


    template<typename DataType>
    matrix<DataType>::matrix(const_matrix_ref m) :
      _stor(0),
      _base(0), _ldim(-1), _nrow(-1), _ncol(-1)
    {
      if (m._base)
        {
          _stor = new data_type[m._nrow * m._ncol];
          _base = _stor;
          _ldim = m._nrow;
          _nrow = m._nrow;
          _ncol = m._ncol;

          copy(m);
        }
    }


    template<typename DataType>
    matrix<DataType>::~matrix()
    {
      if (_stor)
        {
          delete [] _stor;
        }
    }


    template<typename DataType>
    typename matrix<DataType>::matrix_ref
    matrix<DataType>::operator =(data_type x)
    {
      copy(x);

      return *this;
    }


    template<typename DataType>
    typename matrix<DataType>::matrix_ref
    matrix<DataType>::operator =(matrix_type m)
    {
      move(m);

      return *this;
    }


    template<typename DataType>
    typename matrix<DataType>::data_ref
    matrix<DataType>::operator ()(int i, int j)
    {
      MCS_ASSERT((i >= 0) && (i < _nrow));
      MCS_ASSERT((j >= 0) && (j < _ncol));

      return *(_base + i + j * _ldim);
    }


    template<typename DataType>
    typename matrix<DataType>::const_data_ref
    matrix<DataType>::operator ()(int i, int j) const
    {
      MCS_ASSERT((i >= 0) && (i < _nrow));
      MCS_ASSERT((j >= 0) && (j < _ncol));

      return *(_base + i + j * _ldim);
    }


    template<typename DataType>
    typename matrix<DataType>::subvector_ref
    matrix<DataType>::operator ()(int i, range rj)
    {
      MCS_ASSERT((i >= 0) && (i < _nrow));
      MCS_ASSERT((rj._start <= _ncol) && (rj._end <= _ncol));

      if (rj._end == -1)
        {
          rj._end = _ncol;
        }

      return *(new subvector_type(_base + i + rj._start * _ldim,
                                  _ldim,
                                  rj._end - rj._start));
    }


    template<typename DataType>
    typename matrix<DataType>::const_subvector_ref
    matrix<DataType>::operator ()(int i, range rj) const
    {
      MCS_ASSERT((i >= 0) && (i < _nrow));
      MCS_ASSERT((rj._start <= _ncol) && (rj._end <= _ncol));

      if (rj._end == -1)
        {
          rj._end = _ncol;
        }

      return *(new subvector_type(_base + i + rj._start * _ldim,
                                  _ldim,
                                  rj._end - rj._start));
    }


    template<typename DataType>
    typename matrix<DataType>::subvector_ref
    matrix<DataType>::operator ()(range ri, int j)
    {
      MCS_ASSERT((ri._start <= _nrow) && (ri._end <= _nrow));
      MCS_ASSERT((j >= 0) && (j < _ncol));

      if (ri._end == -1)
        {
          ri._end = _nrow;
        }

      return *(new subvector_type(_base + ri._start + j * _ldim,
                                  1,
                                  ri._end - ri._start));
    }


    template<typename DataType>
    typename matrix<DataType>::const_subvector_ref
    matrix<DataType>::operator ()(range ri, int j) const
    {
      MCS_ASSERT((ri._start <= _nrow) && (ri._end <= _nrow));
      MCS_ASSERT((j >= 0) && (j < _ncol));

      if (ri._end == -1)
        {
          ri._end = _nrow;
        }

      return *(new subvector_type(_base + ri._start + j * _ldim,
                                  1,
                                  ri._end - ri._start));
    }


    template<typename DataType>
    typename matrix<DataType>::submatrix_ref
    matrix<DataType>::operator ()(range ri, range rj)
    {
      MCS_ASSERT((ri._start <= _nrow) && (ri._end <= _nrow));
      MCS_ASSERT((rj._start <= _ncol) && (rj._end <= _ncol));

      if (ri._end == -1)
        {
          ri._end = _nrow;
        }

      if (rj._end == -1)
        {
          rj._end = _ncol;
        }

      return *(new submatrix_type(_base + ri._start + rj._start * _ldim,
                                  _ldim,
                                  ri._end - ri._start,
                                  rj._end - rj._start));
    }


    template<typename DataType>
    typename matrix<DataType>::const_submatrix_ref
    matrix<DataType>::operator ()(range ri, range rj) const
    {
      MCS_ASSERT((ri._start <= _nrow) && (ri._end <= _nrow));
      MCS_ASSERT((rj._start <= _ncol) && (rj._end <= _ncol));

      if (ri._end == -1)
        {
          ri._end = _nrow;
        }

      if (rj._end == -1)
        {
          rj._end = _ncol;
        }

      return *(new submatrix_type(_base + ri._start + rj._start * _ldim,
                                  _ldim,
                                  ri._end - ri._start,
                                  rj._end - rj._start));
    }


    template<typename DataType>
    int
    matrix<DataType>::ldim() const
    {
      return _ldim;
    }


    template<typename DataType>
    int
    matrix<DataType>::nrow() const
    {
      return _nrow;
    }


    template<typename DataType>
    int
    matrix<DataType>::ncol() const
    {
      return _ncol;
    }


    template<typename DataType>
    typename matrix<DataType>::subvector_ref
    matrix<DataType>::row(int i)
    {
      MCS_ASSERT((i >= 0) && (i < _nrow));

      return *(new subvector_type(_base + i,
                                  _ldim,
                                  _ncol));
    }


    template<typename DataType>
    typename matrix<DataType>::const_subvector_ref
    matrix<DataType>::row(int i) const
    {
      MCS_ASSERT((i >= 0) && (i < _nrow));

      return *(new subvector_type(_base + i,
                                  _ldim,
                                  _ncol));
    }


    template<typename DataType>
    typename matrix<DataType>::submatrix_ref
    matrix<DataType>::row(range ri)
    {
      MCS_ASSERT((ri._start <= _nrow) && (ri._end <= _nrow));

      if (ri._end == -1)
        {
          ri._end = _nrow;
        }

      return *(new submatrix_type(_base + ri._start,
                                  _ldim,
                                  ri._end - ri._start,
                                  _ncol));
    }


    template<typename DataType>
    typename matrix<DataType>::const_submatrix_ref
    matrix<DataType>::row(range ri) const
    {
      MCS_ASSERT((ri._start <= _nrow) && (ri._end <= _nrow));

      if (ri._end == -1)
        {
          ri._end = _nrow;
        }

      return *(new submatrix_type(_base + ri._start,
                                  _ldim,
                                  ri._end - ri._start,
                                  _ncol));
    }


    template<typename DataType>
    typename matrix<DataType>::subvector_ref
    matrix<DataType>::col(int j)
    {
      MCS_ASSERT((j >= 0) && (j < _ncol));

      return *(new subvector_type(_base + j * _ldim,
                                  1,
                                  _nrow));
    }


    template<typename DataType>
    typename matrix<DataType>::const_subvector_ref
    matrix<DataType>::col(int j) const
    {
      MCS_ASSERT((j >= 0) && (j < _ncol));

      return *(new subvector_type(_base + j * _ldim,
                                  1,
                                  _nrow));
    }


    template<typename DataType>
    typename matrix<DataType>::submatrix_ref
    matrix<DataType>::col(range rj)
    {
      MCS_ASSERT((rj._start <= _ncol) && (rj._end <= _ncol));

      if (rj._end == -1)
        {
          rj._end = _ncol;
        }

      return *(new submatrix_type(_base + rj._start * _ldim,
                                  _ldim,
                                  _nrow,
                                  rj._end - rj._start));
    }


    template<typename DataType>
    typename matrix<DataType>::const_submatrix_ref
    matrix<DataType>::col(range rj) const
    {
      MCS_ASSERT((rj._start <= _ncol) && (rj._end <= _ncol));

      if (rj._end == -1)
        {
          rj._end = _ncol;
        }

      return *(new submatrix_type(_base + rj._start * _ldim,
                                  _ldim,
                                  _nrow,
                                  rj._end - rj._start));
    }


    template<typename DataType>
    typename matrix<DataType>::subvector_ref
    matrix<DataType>::diag()
    {
      return *(new subvector_type(_base,
                                  _ldim + 1,
                                  std::min(_nrow, _ncol)));
    }


    template<typename DataType>
    typename matrix<DataType>::const_subvector_ref
    matrix<DataType>::diag() const
    {
      return *(new subvector_type(_base,
                                  _ldim + 1,
                                  std::min(_nrow, _ncol)));
    }


    template<typename DataType>
    void
    matrix<DataType>::swap(matrix_ref m)
    {
      move(m);
    }


    template<typename DataType>
    void
    matrix<DataType>::swap(submatrix_ref m)
    {
      MCS_ASSERT(_nrow == m._nrow);
      MCS_ASSERT(_ncol == m._ncol);

      data_ptr ptr1 = _base;
      data_ptr ptr2 = m._base;
      for (int j = 0; j < _ncol; ++j)
        {
          for (int i = 0; i < _nrow; ++i)
            {
              std::swap(ptr1[i], ptr2[i]);
            }
          ptr1 += _ldim;
          ptr2 += m._ldim;
        }
    }


    template<typename DataType>
    void
    matrix<DataType>::copy(data_type x)
    {
      data_ptr dst = _base;
      for (int j = 0; j < _ncol; ++j)
        {
          for (int i = 0; i < _nrow; ++i)
            {
              dst[i] = x;
            }
          dst += _ldim;
        }
    }


    template<typename DataType>
    void
    matrix<DataType>::copy(const_matrix_ref m)
    {
      data_ptr src = m._base;
      data_ptr dst = _base;
      for (int j = 0; j < _ncol; ++j)
        {
          for (int i = 0; i < _nrow; ++i)
            {
              dst[i] = src[i];
            }
          src += m._ldim;
          dst += _ldim;
        }
    }


    template<typename DataType>
    void
    matrix<DataType>::move(matrix_ref m)
    {
      std::swap(_stor, m._stor);
      std::swap(_base, m._base);
      std::swap(_ldim, m._ldim);
      std::swap(_nrow, m._nrow);
      std::swap(_ncol, m._ncol);
    }


  }

}


#endif
