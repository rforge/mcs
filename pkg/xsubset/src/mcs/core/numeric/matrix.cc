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
      MATRIX::matrix() :
	alloc_(),
	stor_(0),
	base_(0),
	ldim_(0),
	nrow_(0),
	ncol_(0)
      {
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix(const size_type nrow,
		     const size_type ncol) :
	alloc_(),
	stor_(alloc_.allocate(nrow * ncol)),
	base_(stor_),
	ldim_(nrow),
	nrow_(nrow),
	ncol_(ncol)
      {
	MCS_ASSERT(nrow >= 0);
	MCS_ASSERT(ncol >= 0);
     }


      template<typename Value,
	       typename Alloc>
      template<typename V>
      MATRIX::matrix(const size_type nrow,
		     const size_type ncol,
		     const V val) :
	alloc_(),
	stor_(alloc_.allocate(nrow * ncol)),
	base_(stor_),
	ldim_(nrow),
	nrow_(nrow),
	ncol_(ncol)
      {
	MCS_ASSERT(nrow >= 0);
	MCS_ASSERT(ncol >= 0);

	fill(val);
      }


      template<typename Value,
	       typename Alloc>
      template<typename V>
      MATRIX::matrix(const size_type nrow,
		     const size_type ncol,
		     const V* data) :
	alloc_(),
	stor_(alloc_.allocate(nrow * ncol)),
	base_(stor_),
	ldim_(nrow),
	nrow_(nrow),
	ncol_(ncol)
      {
	MCS_ASSERT(nrow >= 0);
	MCS_ASSERT(ncol >= 0);

	fill(data);
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix(const pointer stor,
		     const size_type nrow,
		     const size_type ncol) :
	alloc_(),
	stor_(stor),
	base_(stor),
	ldim_(nrow),
	nrow_(nrow),
	ncol_(ncol)
      {
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix(const pointer base,
		     const size_type ldim,
		     const size_type nrow,
		     const size_type ncol) :
	alloc_(),
	stor_(0),
	base_(base),
	ldim_(ldim),
	nrow_(nrow),
	ncol_(ncol)
      {
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix(const matrix& mat) :
	alloc_(),
	stor_(alloc_.allocate(mat.nrow_ * mat.ncol_)),
	base_(stor_),
	ldim_(mat.nrow_),
	nrow_(mat.nrow_),
	ncol_(mat.ncol_)
      {
	copy(mat);
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix(matrix&& mat) :
	alloc_(),
	stor_(0),
	base_(0),
	ldim_(1),
	nrow_(0),
	ncol_(0)
      {
	swap(mat);
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::matrix(submatrix_type&& mat) :
	alloc_(),
	stor_(alloc_.allocate(mat.nrow_ * mat.ncol_)),
	base_(stor_),
	ldim_(mat.nrow_),
	nrow_(mat.nrow_),
	ncol_(mat.ncol_)
      {
	copy(mat);
      }


      template<typename Value,
	       typename Alloc>
      MATRIX::~matrix()
      {
	if (stor_)
	  {
	    alloc_.deallocate(stor_, nrow_ * ncol_);
	  }
      }


      template<typename Value,
	       typename Alloc>
      MATRIX&
      MATRIX::operator =(const matrix& mat)
      {
	MCS_ASSERT(mat.nrow_ == nrow_);
	MCS_ASSERT(mat.ncol_ == ncol_);

	copy(mat);

	return *this;
      }


      template<typename Value,
	       typename Alloc>
      MATRIX&
      MATRIX::operator =(matrix&& mat)
      {
	swap(mat);

	return *this;
      }


      template<typename Value,
	       typename Alloc>
      MATRIX&
      MATRIX::operator =(submatrix_type&& mat)
      {
	MCS_ASSERT(mat.nrow_ == nrow_);
	MCS_ASSERT(mat.ncol_ == ncol_);

	copy(mat);

	return *this;
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::reference
      MATRIX::operator ()(const size_type irow,
			  const size_type icol)
      {
	MCS_ASSERT((irow >= 0) && (irow < nrow_));
	MCS_ASSERT((icol >= 0) && (icol < ncol_));

	return at(irow, icol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::const_reference
      MATRIX::operator ()(const size_type irow,
			  const size_type icol) const
      {
	MCS_ASSERT((irow >= 0) && (irow < nrow_));
	MCS_ASSERT((icol >= 0) && (icol < ncol_));

	return at(irow, icol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::subvector_type
      MATRIX::operator ()(const size_type irow,
			  const range_type rcol)
      {
	const size_type icol = rcol.pos();
	const size_type ncol = rcol.open()? ncol_ - icol: rcol.len();

	MCS_ASSERT((irow >= 0) && (irow < nrow_));
	MCS_ASSERT((icol + ncol) <= ncol_);

	return subvector_type(&at(irow, icol), ldim_, ncol);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::subvector_type
      MATRIX::operator ()(const size_type irow,
			  const range_type rcol) const
      {
	const size_type icol = rcol.pos();
	const size_type ncol = rcol.open()? ncol_ - icol: rcol.len();

	MCS_ASSERT((irow >= 0) && (irow < nrow_));
	MCS_ASSERT((icol + ncol) <= ncol_);

	return subvector_type(&at(irow, icol), ldim_, ncol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::subvector_type
      MATRIX::operator ()(const range_type rrow,
			  const size_type icol)
      {
	const size_type irow = rrow.pos();
	const size_type nrow = rrow.open()? nrow_ - irow: rrow.len();

	MCS_ASSERT((irow + nrow) <= nrow_);
	MCS_ASSERT((icol >= 0) && (icol < ncol_));

	return subvector_type(&at(irow, icol), 1, nrow);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::subvector_type
      MATRIX::operator ()(const range_type rrow,
			  const size_type icol) const
      {
	const size_type irow = rrow.pos();
	const size_type nrow = rrow.open()? nrow_ - irow: rrow.len();

	MCS_ASSERT((irow + nrow) <= nrow_);
	MCS_ASSERT((icol >= 0) && (icol < ncol_));

	return subvector_type(&at(irow, icol), 1, nrow);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::submatrix_type
      MATRIX::operator ()(const range_type rrow,
			  const range_type rcol)
      {
	const size_type irow = rrow.pos();
	const size_type nrow = rrow.open()? nrow_ - irow: rrow.len();

	const size_type icol = rcol.pos();
	const size_type ncol = rcol.open()? ncol_ - icol: rcol.len();

	MCS_ASSERT((irow + nrow) <= nrow_);
	MCS_ASSERT((icol + ncol) <= ncol_);

	return submatrix_type(&at(irow, icol), ldim_, nrow, ncol);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::submatrix_type
      MATRIX::operator ()(const range_type rrow,
			  const range_type rcol) const
      {
	const size_type irow = rrow.pos();
	const size_type nrow = rrow.open()? nrow_ - irow: rrow.len();

	const size_type icol = rcol.pos();
	const size_type ncol = rcol.open()? ncol_ - icol: rcol.len();

	MCS_ASSERT((irow + nrow) <= nrow_);
	MCS_ASSERT((icol + ncol) <= ncol_);

	return submatrix_type(&at(irow, icol), ldim_, nrow, ncol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::size_type
      MATRIX::ldim() const
      {
	return ldim_;
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::size_type
      MATRIX::nrow() const
      {
	return nrow_;
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::size_type
      MATRIX::ncol() const
      {
	return ncol_;
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::subvector_type
      MATRIX::row(const size_type irow)
      {
	MCS_ASSERT((irow >= 0) && (irow < nrow_));

	return subvector_type(&at(irow, 0), ldim_, ncol_);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::subvector_type
      MATRIX::row(const size_type irow) const
      {
	MCS_ASSERT((irow >= 0) && (irow < nrow_));

	return subvector_type(&at(irow, 0), ldim_, ncol_);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::subvector_type
      MATRIX::col(const size_type icol)
      {
	MCS_ASSERT((icol >= 0) && (icol < ncol_));

	return subvector_type(&at(0, icol), 1, nrow_);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::subvector_type
      MATRIX::col(size_type icol) const
      {
	MCS_ASSERT((icol >= 0) && (icol < ncol_));

	return subvector_type(&at(0, icol), 1, nrow_);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::submatrix_type
      MATRIX::row(const range_type rrow)
      {
	const size_type irow = rrow.pos();
	const size_type nrow = rrow.open()? nrow_ - irow: rrow.len();

	MCS_ASSERT((irow + nrow) <= nrow_);

	return submatrix_type(&at(irow, 0), ldim_, nrow, ncol_);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::submatrix_type
      MATRIX::row(const range_type rrow) const
      {
	const size_type irow = rrow.pos();
	const size_type nrow = rrow.open()? nrow_ - irow: rrow.len();

	MCS_ASSERT((irow + nrow) <= nrow_);

	return submatrix_type(&at(irow, 0), ldim_, nrow, ncol_);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::submatrix_type
      MATRIX::col(const range_type rcol)
      {
	const size_type icol = rcol.pos();
	const size_type ncol = rcol.open()? ncol_ - icol: rcol.len();

	MCS_ASSERT((icol + ncol) <= ncol_);

	return submatrix_type(&at(0, icol), ldim_, nrow_, ncol);
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::submatrix_type
      MATRIX::col(const range_type rcol) const
      {
	const size_type icol = rcol.pos();
	const size_type ncol = rcol.open()? ncol_ - icol: rcol.len();

	MCS_ASSERT((icol + ncol) <= ncol_);

	return submatrix_type(at(0, icol), ldim_, nrow_, ncol);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::subvector_type
      MATRIX::diag()
      {
	return subvector_type(base_, ldim_ + 1, std::min(nrow_, ncol_));
      }


      template<typename Value,
	       typename Alloc>
      const typename MATRIX::subvector_type
      MATRIX::diag() const
      {
	return subvector_type(base_, ldim_ + 1, std::min(nrow_, ncol_));
      }


      template<typename Value,
	       typename Alloc>
      template<typename V>
      void
      MATRIX::fill(const V val)
      {
	pointer ptr = base_;
	for (size_type j = 0; j < ncol_; ++j)
	  {
	    for (size_type i = 0; i < nrow_; ++i)
	      {
		ptr[i] = static_cast<value_type>(val);
	      }
	    ptr += ldim_;
	  }
      }


      template<typename Value,
	       typename Alloc>
      template<typename V>
      void
      MATRIX::fill(const V* data)
      {
	pointer ptr = base_;
	for (size_type j = 0; j < ncol_; ++j)
	  {
	    for (size_type i = 0; i < nrow_; ++i)
	      {
		ptr[i] = static_cast<value_type>(*data++);
	      }
	    ptr += ldim_;
	  }
      }


      template<typename Value,
	       typename Alloc>
      void
      MATRIX::copy(const matrix& mat)
      {
	MCS_ASSERT(nrow_ == mat.nrow_);
	MCS_ASSERT(ncol_ == mat.ncol_);

	const size_type lsrc = mat.ldim_;
	const size_type ldst = ldim_;
	const size_type nrow = nrow_;
	const size_type ncol = ncol_;

	const_pointer src = mat.base_;
	pointer dst = base_;
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
      MATRIX::swap(matrix& mat)
      {
        // TODO: swap allocator
	std::swap(stor_, mat.stor_);
	std::swap(base_, mat.base_);
	std::swap(ldim_, mat.ldim_);
	std::swap(nrow_, mat.nrow_);
	std::swap(ncol_, mat.ncol_);
      }


      template<typename Value,
	       typename Alloc>
      typename MATRIX::reference
      MATRIX::at(const size_type irow,
		 const size_type icol) const
      {
	return base_[irow + icol * ldim_];
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
      SUBMATRIX::operator =(const matrix_type& mat)
      {
	MCS_ASSERT(mat.nrow_ == nrow_);
	MCS_ASSERT(mat.ncol_ == ncol_);

	copy(mat);

	return *this;
      }


      template<typename Value,
	       typename Alloc>
      SUBMATRIX&
      SUBMATRIX::operator =(matrix_type&& mat)
      {
	MCS_ASSERT(mat.nrow_ == nrow_);
	MCS_ASSERT(mat.ncol_ == ncol_);

	copy(mat);

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
