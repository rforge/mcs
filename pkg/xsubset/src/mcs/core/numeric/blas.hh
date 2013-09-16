/**
 * @file blas.hh
 */
#ifndef MCS_CORE_NUMERIC_BLAS_HH
#define MCS_CORE_NUMERIC_BLAS_HH


#include "traits.hh"


namespace mcs     {
namespace core    {
namespace numeric {


template<typename Value>
struct blas
{
  typedef typename traits<Value>::size_type size_type;
  typedef typename traits<Value>::value_type value_type;
  typedef typename traits<Value>::reference reference;
  typedef typename traits<Value>::const_reference const_reference;
  typedef typename traits<Value>::pointer pointer;
  typedef typename traits<Value>::const_pointer const_pointer;

  typedef typename traits<Value>::vector vector;
  typedef typename traits<Value>::const_vector const_vector;
  typedef typename traits<Value>::vector_reference vector_reference;
  typedef typename traits<Value>::const_vector_reference const_vector_reference;

  typedef typename traits<Value>::matrix matrix;
  typedef typename traits<Value>::const_matrix const_matrix;
  typedef typename traits<Value>::matrix_reference matrix_reference;
  typedef typename traits<Value>::const_matrix_reference const_matrix_reference;


  void copy(size_type n, const_reference x, size_type incx,
	    reference y, size_type incy);

  void copy(const_vector_reference_type x, vector_reference_type y)
  {
    copy(x.len(), x(0), x.inc(), y(0), y.inc());
  }

  void copy(subscript s, const_vector_reference x, vector_reference y)
  {
    copy(x(s), y(s));
  }

  void gemv(const char* trans, size_type m, size_type n,
	    value_type alpha, const_reference a, size_type lda,
	    const_reference x, size_type incx, value_type beta,
	    reference y, size_type incy);

  void gemv(const char* trans, value_type alpha, const_matrix_reference_type a,
	    const_vector_reference_type x, value_type beta,
            vector_reference_type y)
  {
    gemv(trans, a.nrow(), a.ncol(), alpha, a(0, 0), a.ldim(),
	 x(0), x.inc(), beta, y(0), y.inc());
  }

  void rotg(reference x, reference y, reference c, reference s);

  void rot(size_type n, reference x, size_type incx, reference y, size_type incy,
	   value_type c, value_type s);

  void rot(vector_reference_type x, vector_reference_type y,
           value_type c, value_type s)
  {
    rot(x.len(), x(0), x.inc(), y(0), y.inc(), c, s);
  }

  void rot(vector_reference_type x, vector_reference_type y)
  {
    value_type c, s;

    rotg(x(0), y(0), c, s);
    rot(x, y, c, s);
  }

  void rot(subscript s, vector_reference x, vector_reference y)
  {
    rot(x(s), y(s));
  }


}  // numeric
}  // core
}  // mcs


#include "blas.cc"
#endif
