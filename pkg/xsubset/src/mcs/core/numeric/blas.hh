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
  typedef typename traits<Value>::reference_type reference_type;
  typedef typename traits<Value>::const_reference_type const_reference_type;
  typedef typename traits<Value>::pointer_type pointer_type;
  typedef typename traits<Value>::const_pointer_type const_pointer_type;

  typedef typename traits<Value>::vector_type vector_type;
  typedef typename traits<Value>::const_vector_type const_vector_type;
  typedef typename traits<Value>::vector_reference_type vector_reference_type;
  typedef typename traits<Value>::const_vector_reference_type const_vector_reference_type;

  typedef typename traits<Value>::matrix_type matrix_type;
  typedef typename traits<Value>::const_matrix_type const_matrix_type;
  typedef typename traits<Value>::matrix_reference_type matrix_reference_type;
  typedef typename traits<Value>::const_matrix_reference_type const_matrix_reference_type;


  void copy(size_type n, const_reference_type x, size_type incx,
	    reference_type y, size_type incy);

  void copy(const_vector_reference_type x, vector_reference_type y)
  {
    copy(x.len(), x(0), x.inc(), y(0), y.inc());
  }

  void copy(subscript_type s, const_vector_reference_type x, vector_reference_type y)
  {
    copy(x(s), y(s));
  }

  void gemv(const char* trans, size_type m, size_type n,
	    value_type alpha, const_reference_type a, size_type lda,
	    const_reference_type x, size_type incx, value_type beta,
	    reference_type y, size_type incy);

  void gemv(const char* trans, value_type alpha, const_matrix_reference_type a,
	    const_vector_reference_type x, value_type beta,
            vector_reference_type y)
  {
    gemv(trans, a.nrow(), a.ncol(), alpha, a(0, 0), a.ldim(),
	 x(0), x.inc(), beta, y(0), y.inc());
  }

  void rotg(reference_type x, reference_type y,
            reference_type c, reference_type s);

  void rot(size_type n, reference_type x, size_type incx, reference_type y,
           size_type incy, value_type c, value_type s);

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

  void rot(subscript_type s, vector_reference_type x, vector_reference_type y)
  {
    rot(x(s), y(s));
  }


}  // numeric
}  // core
}  // mcs


#include "blas.cc"
#endif
