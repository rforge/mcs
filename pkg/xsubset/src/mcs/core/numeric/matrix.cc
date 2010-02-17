/**
 * @file matrix.cc
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_NUMERIC_MATRIX_CC
#define MCS_CORE_NUMERIC_MATRIX_CC


#include <algorithm>

#include "../../mcs.hh"

#include "range.hh"
#include "vector.hh"
#include "matrix.hh"


#define MATRIX matrix<Value, Alloc>
#define SUBMATRIX submatrix<Value, Alloc>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix_impl::matrix_impl() :
	allocator_type(),
	stor(0),
	base(0),
	ldim(0),
	nrow(0),
	ncol(0)
      {
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix_impl::matrix_impl(const allocator_type& alloc):
	allocator_type(alloc),
	stor(0),
	base(0),
	ldim(0),
	nrow(0),
	ncol(0)
      {
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix() :
	impl_()
      {
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix(const size_type nrow,
		     const size_type ncol) :
	impl_()
      {
	MCS_ASSERT(nrow >= 0);
	MCS_ASSERT(ncol >= 0);

	impl_.stor = impl_.allocate(nrow * ncol);
	impl_.base = impl_.stor;
	impl_.ldim = nrow;
	impl_.nrow = nrow;
	impl_.ncol = ncol;
      }


      template<typename Value,
	       typename Alloc>
      template<typename V>
      MATRIX::matrix(const size_type nrow,
		     const size_type ncol,
		     const V val) :
	impl_()
      {
	MCS_ASSERT(nrow >= 0);
	MCS_ASSERT(ncol >= 0);

	impl_.stor = impl_.allocate(nrow * ncol);
	impl_.base = impl_.stor;
	impl_.ldim = nrow;
	impl_.nrow = nrow;
	impl_.ncol = ncol;

	fill(val);
      }


      template<typename Value,
	       typename Alloc>
      template<typename V>
      MATRIX::matrix(const size_type nrow,
		     const size_type ncol,
		     const V* data) :
	impl_()
      {
	MCS_ASSERT(nrow >= 0);
	MCS_ASSERT(ncol >= 0);

	const size_type size = nrow * ncol;

	impl_.stor = impl_.allocate(size);
	impl_.base = impl_.stor;
	impl_.ldim = nrow;
	impl_.nrow = nrow;
	impl_.ncol = ncol;

	fill(data);
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix(const pointer stor,
		     const size_type nrow,
		     const size_type ncol) :
	impl_()
      {
	impl_.stor = stor;
	impl_.base = stor;
	impl_.ldim = nrow;
	impl_.nrow = nrow;
	impl_.ncol = ncol;
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix(const pointer base,
		     const size_type ldim,
		     const size_type nrow,
		     const size_type ncol) :
	impl_()
      {
	impl_.stor = 0;
	impl_.base = base;
	impl_.ldim = ldim;
	impl_.nrow = nrow;
	impl_.ncol = ncol;
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix(const MATRIX& mat) :
	impl_()
      {
	if (mat.impl_.base)
	  {
	    impl_.stor = impl_.allocate(mat.impl_.nrow * mat.impl_.ncol);
	    impl_.base = impl_.stor;
	    impl_.ldim = mat.impl_.nrow;
	    impl_.nrow = mat.impl_.nrow;
	    impl_.ncol = mat.impl_.ncol;

	    copy(mat);
	  }
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix(MATRIX&& mat) :
	impl_()
      {
	swap(mat);
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix(submatrix_type&& mat) :
	impl_()
      {
	impl_.stor = impl_.allocate(mat.impl_.nrow * mat.impl_.ncol);
	impl_.base = impl_.stor;
	impl_.ldim = mat.impl_.nrow;
	impl_.nrow = mat.impl_.nrow;
	impl_.ncol = mat.impl_.ncol;

	copy(mat);
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::~matrix()
      {
	if (impl_.stor)
	  {
	    impl_.deallocate(impl_.stor, 0);
	  }
      }


      template<typename Value,
	       typename Alloc>
      MATRIX&
      MATRIX::operator =(MATRIX mat)
      {
	swap(mat);

	return *this;
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::reference
      MATRIX::operator ()(const size_type irow,
			  const size_type icol)
      {
	MCS_ASSERT((irow >= 0) && (irow < impl_.nrow));
	MCS_ASSERT((icol >= 0) && (icol < impl_.ncol));

	return at(irow, icol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::const_reference
      MATRIX::operator ()(const size_type irow,
			  const size_type icol) const
      {
	MCS_ASSERT((irow >= 0) && (irow < impl_.nrow));
	MCS_ASSERT((icol >= 0) && (icol < impl_.ncol));

	return at(irow, icol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::subvector_type
      MATRIX::operator ()(const size_type irow,
			  const range_type rcol)
      {
	const size_type icol = rcol.pos();
	const size_type ncol = rcol.open()? impl_.ncol - icol: rcol.len();

	MCS_ASSERT((irow >= 0) && (irow < impl_.nrow));
	MCS_ASSERT((icol + ncol) <= impl_.ncol);

	return subvector_type(&at(irow, icol), impl_.ldim, ncol);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::subvector_type
      MATRIX::operator ()(const size_type irow,
			  const range_type rcol) const
      {
	const size_type icol = rcol.pos();
	const size_type ncol = rcol.open()? impl_.ncol - icol: rcol.len();

	MCS_ASSERT((irow >= 0) && (irow < impl_.nrow));
	MCS_ASSERT((icol + ncol) <= impl_.ncol);

	return subvector_type(&at(irow, icol), impl_.ldim, ncol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::subvector_type
      MATRIX::operator ()(const range_type rrow,
			  const size_type icol)
      {
	const size_type irow = rrow.pos();
	const size_type nrow = rrow.open()? impl_.nrow - irow: rrow.len();

	MCS_ASSERT((irow + nrow) <= impl_.nrow);
	MCS_ASSERT((icol >= 0) && (icol < impl_.ncol));

	return subvector_type(&at(irow, icol), 1, nrow);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::subvector_type
      MATRIX::operator ()(const range_type rrow,
			  const size_type icol) const
      {
	const size_type irow = rrow.pos();
	const size_type nrow = rrow.open()? impl_.nrow - irow: rrow.len();

	MCS_ASSERT((irow + nrow) <= impl_.nrow);
	MCS_ASSERT((icol >= 0) && (icol < impl_.ncol));

	return subvector_type(&at(irow, icol), 1, nrow);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::submatrix_type
      MATRIX::operator ()(const range_type rrow,
			  const range_type rcol)
      {
	const size_type irow = rrow.pos();
	const size_type nrow = rrow.open()? impl_.nrow - irow: rrow.len();

	const size_type icol = rcol.pos();
	const size_type ncol = rcol.open()? impl_.ncol - icol: rcol.len();

	MCS_ASSERT((irow + nrow) <= impl_.nrow);
	MCS_ASSERT((icol + ncol) <= impl_.ncol);

	return submatrix_type(&at(irow, icol), impl_.ldim, nrow, ncol);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::submatrix_type
      MATRIX::operator ()(const range_type rrow,
			  const range_type rcol) const
      {
	const size_type irow = rrow.pos();
	const size_type nrow = rrow.open()? impl_.nrow - irow: rrow.len();

	const size_type icol = rcol.pos();
	const size_type ncol = rcol.open()? impl_.ncol - icol: rcol.len();

	MCS_ASSERT((irow + nrow) <= impl_.nrow);
	MCS_ASSERT((icol + ncol) <= impl_.ncol);

	return submatrix_type(&at(irow, icol), impl_.ldim, nrow, ncol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::size_type
      MATRIX::ldim() const
      {
	return impl_.ldim;
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::size_type
      MATRIX::nrow() const
      {
	return impl_.nrow;
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::size_type
      MATRIX::ncol() const
      {
	return impl_.ncol;
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::subvector_type
      MATRIX::row(const size_type irow)
      {
	MCS_ASSERT((irow >= 0) && (irow < impl_.nrow));

	return subvector_type(&at(irow, 0), impl_.ldim, impl_.ncol);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::subvector_type
      MATRIX::row(const size_type irow) const
      {
	MCS_ASSERT((irow >= 0) && (irow < impl_.nrow));

	return subvector_type(&at(irow, 0), impl_.ldim, impl_.ncol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::subvector_type
      MATRIX::col(const size_type icol)
      {
	MCS_ASSERT((icol >= 0) && (icol < impl_.ncol));

	return subvector_type(&at(0, icol), 1, impl_.nrow);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::subvector_type
      MATRIX::col(size_type icol) const
      {
	MCS_ASSERT((icol >= 0) && (icol < impl_.ncol));

	return subvector_type(&at(0, icol), 1, impl_.nrow);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::submatrix_type
      MATRIX::row(const range_type rrow)
      {
	const size_type irow = rrow.pos();
	const size_type nrow = rrow.open()? impl_.nrow - irow: rrow.len();

	MCS_ASSERT((irow + nrow) <= impl_.nrow);

	return submatrix_type(&at(irow, 0), impl_.ldim, nrow, impl_.ncol);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::submatrix_type
      MATRIX::row(const range_type rrow) const
      {
	const size_type irow = rrow.pos();
	const size_type nrow = rrow.open()? impl_.nrow - irow: rrow.len();

	MCS_ASSERT((irow + nrow) <= impl_.nrow);

	return submatrix_type(&at(irow, 0), impl_.ldim, nrow, impl_.ncol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::submatrix_type
      MATRIX::col(const range_type rcol)
      {
	const size_type icol = rcol.pos();
	const size_type ncol = rcol.open()? impl_.ncol - icol: rcol.len();

	MCS_ASSERT((icol + ncol) <= impl_.ncol);

	return submatrix_type(&at(0, icol), impl_.ldim, impl_.nrow, ncol);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::submatrix_type
      MATRIX::col(const range_type rcol) const
      {
	const size_type icol = rcol.pos();
	const size_type ncol = rcol.open()? impl_.ncol - icol: rcol.len();

	MCS_ASSERT((icol + ncol) <= impl_.ncol);

	return submatrix_type(at(0, icol), impl_.ldim, impl_.nrow, ncol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::subvector_type
      MATRIX::diag()
      {
	return subvector_type(impl_.base, impl_.ldim + 1, std::min(impl_.nrow, impl_.ncol));
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::subvector_type
      MATRIX::diag() const
      {
	return subvector_type(impl_.base, impl_.ldim + 1, std::min(impl_.nrow, impl_.ncol));
      }


      template<typename Value,
	       typename Alloc>
      template<typename V>
      void
      MATRIX::fill(const V val)
      {
	const size_type ldim = impl_.ldim;
	const size_type nrow = impl_.nrow;
	const size_type ncol = impl_.ncol;

	const value_type fill_val = static_cast<value_type>(val);

	pointer dst = impl_.base;
	for (size_type j = 0; j < ncol; ++j)
	  {
	    for (size_type i = 0; i < nrow; ++i)
	      {
		dst[i] = fill_val;
	      }
	    dst += ldim;
	  }
      }


      template<typename Value,
	       typename Alloc>
      template<typename V>
      void
      MATRIX::fill(const V* data)
      {
	const size_type ldim = impl_.ldim;
	const size_type nrow = impl_.nrow;
	const size_type ncol = impl_.ncol;

	pointer dst = impl_.base;
	for (size_type j = 0; j < ncol; ++j)
	  {
	    for (size_type i = 0; i < nrow; ++i)
	      {
		dst[i] = static_cast<value_type>(*data++);
	      }
	    dst += ldim;
	  }
      }


      template<typename Value,
	       typename Alloc>
      void
      MATRIX::copy(const MATRIX& mat)
      {
	MCS_ASSERT(impl_.nrow == mat.impl_.nrow);
	MCS_ASSERT(impl_.ncol == mat.impl_.ncol);

	const size_type lsrc = mat.impl_.ldim;
	const size_type ldst = impl_.ldim;
	const size_type nrow = impl_.nrow;
	const size_type ncol = impl_.ncol;

	const_pointer src = mat.impl_.base;
	pointer dst = impl_.base;
	for (size_type j = 0; j < ncol; ++j)
	  {
	    for (size_type i = 0; i < nrow; ++i)
	      {
		dst[i] = src[i];
	      }
	    src += lsrc;
	    dst += ldst;
	  }
      }


      template<typename Value,
	       typename Alloc>
      void
      MATRIX::swap(MATRIX& mat)
      {
	std::swap(impl_.stor, mat.impl_.stor);
	std::swap(impl_.base, mat.impl_.base);
	std::swap(impl_.ldim, mat.impl_.ldim);
	std::swap(impl_.nrow, mat.impl_.nrow);
	std::swap(impl_.ncol, mat.impl_.ncol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::reference
      MATRIX::at(const size_type irow,
		 const size_type icol) const
      {
	return impl_.base[irow + icol * impl_.ldim];
      }


      template<typename Value,
	       typename Alloc>
      SUBMATRIX::submatrix(const pointer base,
			   const size_type ldim,
			   const size_type nrow,
			   const size_type ncol)
	: matrix_type(base, ldim, nrow, ncol)
      {
      }


      template<typename Value,
	       typename Alloc>
      SUBMATRIX::~submatrix()
      {
      }


      template<typename Value,
	       typename Alloc>
      SUBMATRIX&
      SUBMATRIX::operator =(const MATRIX& mat)
      {
	MCS_ASSERT(mat.nrow() == this->impl_.nrow);
	MCS_ASSERT(mat.ncol() == this->impl_.ncol);

	matrix_type::copy(mat);

	return *this;
      }


      template<typename Value,
	       typename Alloc>
      SUBMATRIX&
      SUBMATRIX::operator =(MATRIX&& mat)
      {
	MCS_ASSERT(mat.nrow() == this->impl_.nrow);
	MCS_ASSERT(mat.ncol() == this->impl_.ncol);

	matrix_type::copy(mat);

	return *this;
      }


    }

  }

}


#undef MATRIX
#undef SUBMATRIX


#define MATRIX mcs::core::numeric::matrix<Value, Alloc>
#define SUBMATRIX mcs::core::numeric::submatrix<Value, Alloc>


namespace std
{


  template<typename Value,
	   typename Alloc>
  void
  swap(MATRIX& x, MATRIX& y)
  {
    x.swap(y);
  }


}


#undef MATRIX
#undef SUBMATRIX


#endif
