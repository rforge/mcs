/**
 * @file lapack.hh
 */
#ifndef MCS_CORE_NUMERIC_LAPACK_HH
#define MCS_CORE_NUMERIC_LAPACK_HH


#include <algorithm>

#include "traits.hh"


namespace mcs     {
namespace core    {
namespace numeric {


template<typename Value>
struct lapack
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

  typedef typename traits<Value>::subscript subscript;


  int info;


  void lacpy(const char* uplo, const size_type m, const size_type n, const_reference a,
	     const size_type lda, reference b, const size_type ldb);

  void lacpy(const char* uplo, const_matrix_reference a, matrix_reference b)
  {
    lacpy(uplo, a.nrow(), a.ncol(), a(0, 0), a.ldim(), b(0, 0), b.ldim());
  }

  void lacpy(const char* uplo, subscript row, subscript col,
	     const_matrix_reference a, matrix_reference b)
  {
    lacpy(uplo, a(row, col), b(row, col));
  }

  void geqrf(const size_type m, const size_type n, reference a, const size_type lda,
	     reference tau, reference work, const size_type lwork);

  void geqrf(matrix_reference a)
  {
    const size_type m = a.nrow();
    const size_type n = a.ncol();
    const size_type lwork = std::max(1, n);
    
    value_type tau[std::min(m, n)];
    value_type work[lwork];

    geqrf(m, n, a(0, 0), a.ldim(), *tau, *work, lwork);
  }

  void geqrf(subscript pos, matrix_reference a)
  {
    geqrf(a(pos, pos));
  }

}


}  // numeric
}  // core
}  // mcs


#include "lapack.cc"
#endif
