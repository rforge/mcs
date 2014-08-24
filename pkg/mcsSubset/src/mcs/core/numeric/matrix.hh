#ifndef _MCS_CORE_NUMERIC_MATRIX_HH_
#define _MCS_CORE_NUMERIC_MATRIX_HH_


#include <algorithm>
#include <cassert>
#include <utility>

#include "../../mcs.hh"

#include "subscript.hh"
#include "traits.hh"
#include "vector.hh"


namespace mcs     {
namespace core    {
namespace numeric {


// matrix base
template<typename Derived>
class matrix_base
{
  template<typename D>
  friend class matrix_base;
public:
  typedef matrix_base<Derived> self_type;
  typedef Derived derived_type;

  typedef typename traits<Derived>::size_type size_type;
  typedef typename traits<Derived>::value_type value_type;
  typedef typename traits<Derived>::reference reference;
  typedef typename traits<Derived>::const_reference const_reference;
  typedef typename traits<Derived>::pointer pointer;
  typedef typename traits<Derived>::const_pointer const_pointer;

  typedef typename traits<Derived>::vector vector;
  typedef typename traits<Derived>::const_vector const_vector;
  typedef typename traits<Derived>::vector_reference vector_reference;
  typedef typename traits<Derived>::const_vector_reference const_vector_reference;

  typedef typename traits<Derived>::matrix matrix;
  typedef typename traits<Derived>::const_matrix const_matrix;
  typedef typename traits<Derived>::matrix_reference matrix_reference;
  typedef typename traits<Derived>::const_matrix_reference const_matrix_reference;

  typedef typename traits<Derived>::subscript subscript;
protected:
  pointer buf_;
  pointer ptr_;
  size_type nrow_;
  size_type ncol_;
  size_type ldim_;
protected:
  matrix_base() :
    buf_(nullptr),
    ptr_(nullptr),
    nrow_(0),
    ncol_(0),
    ldim_(0)
  {
  }
  matrix_base(const size_type nrow, const size_type ncol) :
    buf_(new value_type[nrow * ncol]),
    ptr_(buf_),
    nrow_(nrow),
    ncol_(ncol),
    ldim_(nrow)
  {
    MCS_ASSERT(nrow >= 0, "invalid argument: nrow (matrix_base::matrix_base)");
    MCS_ASSERT(ncol >= 0, "invalid argument: ncol (matrix_base::matrix_base)");
  }
  matrix_base(pointer ptr, const size_type nrow,
	      const size_type ncol, const size_type ldim) :
    buf_(nullptr),
    ptr_(ptr),
    nrow_(nrow),
    ncol_(ncol),
    ldim_(ldim)
  {
    MCS_ASSERT(ptr != nullptr, "invalid argument: ptr (matrix_base::matrix_base)");
    MCS_ASSERT(nrow >= 0, "invalid argument: nrow (matrix_base::matrix_base)");
    MCS_ASSERT(ncol >= 0, "invalid argument: ncol (matrix_base::matrix_base)");
    MCS_ASSERT(ldim >= nrow, "invalid argument: ldim (matrix_base::matrix_base)");
  }
  ~matrix_base()
  {
    if (buf_ != nullptr)
      {
        delete [] buf_;
        buf_ = nullptr;
      }
  }
  reference at(const size_type row, const size_type col)
  {
    MCS_ASSERT((0 <= row) && (row < nrow_), "invalid argument: row (matrix_base::at)");
    MCS_ASSERT((0 <= col) && (col < ncol_), "invalid argument: col (matrix_base::at)");

    return ptr_[col * ldim_ + row];
  }
  const_reference at(const size_type row, const size_type col) const
  {
    MCS_ASSERT((0 <= row) && (row < nrow_), "invalid argument: row (matrix_base::at)");
    MCS_ASSERT((0 <= col) && (col < ncol_), "invalid argument: col (matrix_base::at)");

    return ptr_[col * ldim_ + row];
  }
  matrix_reference at(const subscript& row, const subscript& col)
  {
    MCS_ASSERT(row.inc == 1, "invalid argument: row (matrix_base::copy)");
    MCS_ASSERT(col.inc == 1, "invalid argument: col (matrix_base::copy)");

    return matrix_reference(buf_ + col.pos * ldim_ + row.pos, row.len, col.len, ldim_);
  }
  const_matrix_reference at(const subscript& row, const subscript& col) const
  {
    return const_matrix_reference(buf_ + col.pos * ldim_ + row.pos, row.len, col.len, ldim_);
  }
  vector_reference row(const size_type pos)
  {
    return vector_reference(buf_ + pos, ncol_, ldim_);
  }
  const_vector_reference row(const size_type pos) const
  {
    return const_vector_reference(buf_ + pos, ncol_, ldim_);
  }
  vector_reference col(const size_type pos)
  {
    return vector_reference(buf_ + pos * ldim_, nrow_, 1);
  }
  const_vector_reference col(const size_type pos) const
  {
    return const_vector_reference(buf_ + pos * ldim_, nrow_, 1);
  }
  template<typename D>
  void copy(const matrix_base<D>& x)
  {
    MCS_ASSERT(nrow_ == x.nrow_, "invalid argument: x (matrix_base::copy)");
    MCS_ASSERT(ncol_ == x.ncol_, "invalid argument: x (matrix_base::copy)");

    pointer dst = buf_;
    const_pointer src = x.buf_;

    for (int j = 0; j < ncol_; ++j)
      {
	std::copy(src, src + nrow_, dst);

	src += x.ldim_;
	dst += ldim_;
      }
  }
  template<typename D>
  void move(matrix_base<D>& x)
  {
    std::swap(buf_, x.buf_);
    std::swap(ptr_, x.ptr_);
    std::swap(nrow_, x.nrow_);
    std::swap(ncol_, x.ncol_);
    std::swap(ldim_, x.ldim_);
  }
  template<typename D>
  void alias(const matrix_base<D>& x)
  {
    ptr_ = x.ptr_;
    nrow_ = x.nrow_;
    ncol_ = x.ncol_;
    ldim_ = x.ldim_;
  }
};


// matrix, const_matrix, matrix_reference, const_matrix_reference
template<typename Value>
class matrix : private matrix_base<matrix<Value> >
{
  friend class const_matrix<Value>;
  friend class matrix_reference<Value>;
  friend class const_matrix_reference<Value>;
public:
  typedef matrix<Value> self_type;
  typedef matrix_base<self_type> base_type;

  typedef typename base_type::size_type size_type;
  typedef typename base_type::value_type value;
  typedef typename base_type::reference reference;
  typedef typename base_type::const_reference const_reference;
  typedef typename base_type::pointer pointer;
  typedef typename base_type::const_pointer const_pointer;

  typedef typename base_type::vector vector;
  typedef typename base_type::const_vector const_vector;
  typedef typename base_type::vector_reference vector_reference;
  typedef typename base_type::const_vector_reference const_vector_reference;

  //typedef typename base_type::matrix matrix;
  typedef typename base_type::const_matrix const_matrix;
  typedef typename base_type::matrix_reference matrix_reference;
  typedef typename base_type::const_matrix_reference const_matrix_reference;

  typedef typename base_type::subscript subscript;
public:
  matrix(const size_type nrow, const size_type ncol) :
    base_type(nrow, ncol)
  {
  }
  matrix(const matrix& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  matrix(matrix&& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::move(x);
  }
  matrix(const const_matrix& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  matrix(const_matrix&& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::move(x);
  }
  matrix(const matrix_reference& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  matrix(const const_matrix_reference& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  self_type& operator =(matrix x)
  {
    MCS_ASSERT(base_type::nrow_ == x.nrow_, "invalid argument: x (matrix::operator=)");
    MCS_ASSERT(base_type::ncol_ == x.ncol_, "invalid argument: x (matrix::operator=)");

    base_type::move(x);

    return *this;
  }
  reference operator ()(const size_type row, const size_type col)
  {
    return base_type::at(row, col);
  }
  const_reference operator ()(const size_type row, const size_type col) const
  {
    return base_type::at(row, col);
  }
  matrix_reference operator ()(const subscript& row, const subscript& col)
  {
    return base_type::at(row, col);
  }
  const_matrix_reference operator ()(const subscript& row, const subscript& col) const
  {
    return base_type::at(row, col);
  }
  vector_reference row(const size_type pos)
  {
    base_type::row(pos);
  }
  const_vector_reference row(const size_type pos) const
  {
    base_type::row(pos);
  }
  vector_reference col(const size_type pos)
  {
    base_type::col(pos);
  }
  const_vector_reference col(const size_type pos) const
  {
    base_type::col(pos);
  }
  pointer ptr()
  {
    return base_type::ptr_;
  }
  const_pointer ptr() const
  {
    return base_type::ptr_;
  }
  size_type nrow() const
  {
    return base_type::nrow_;
  }
  size_type ncol() const
  {
    return base_type::ncol_;
  }
  size_type ldim() const
  {
    return base_type::ldim_;
  }
};


template<typename Value>
class const_matrix : private matrix_base<const_matrix<Value> >
{
  friend class matrix<Value>;
  friend class matrix_reference<Value>;
  friend class const_matrix_reference<Value>;
public:
  typedef const_matrix<Value> self_type;
  typedef matrix_base<self_type> base_type;

  typedef typename base_type::size_type size_type;
  typedef typename base_type::value_type value_type;
  typedef typename base_type::reference reference;
  typedef typename base_type::const_reference const_reference;
  typedef typename base_type::pointer pointer;
  typedef typename base_type::const_pointer const_pointer;

  typedef typename base_type::vector vector;
  typedef typename base_type::const_vector const_vector;
  typedef typename base_type::vector_reference vector_reference;
  typedef typename base_type::const_vector_reference const_vector_reference;

  typedef typename base_type::matrix matrix;
  //typedef typename base_type::const_matrix const_matrix;
  typedef typename base_type::matrix_reference matrix_reference;
  typedef typename base_type::const_matrix_reference const_matrix_reference;

  typedef typename base_type::subscript subscript;
public:
  const_matrix(const size_type nrow, const size_type ncol) :
    base_type(nrow, ncol)
  {
  }
  const_matrix(const matrix& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  const_matrix(matrix&& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::move(x);
  }
  const_matrix(const const_matrix& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  const_matrix(const_matrix&& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::move(x);
  }
  const_matrix(const matrix_reference& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  const_matrix(const const_matrix_reference& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  const_reference operator ()(const size_type row, const size_type col) const
  {
    return base_type::at(row, col);
  }
  const_matrix_reference operator ()(const subscript& row, const subscript& col) const
  {
    return base_type::at(row, col);
  }
  const_vector_reference row(const size_type pos) const
  {
    base_type::row(pos);
  }
  const_vector_reference col(const size_type pos) const
  {
    base_type::col(pos);
  }
  const_pointer ptr() const
  {
    return base_type::ptr_;
  }
  size_type nrow() const
  {
    return base_type::rnow_;
  }
  size_type ncol() const
  {
    return base_type::rcol_;
  }
  size_type ldim() const
  {
    return base_type::ldim_;
  }
};


template<typename Value>
class matrix_reference : private matrix_base<matrix_reference<Value> >
{
  friend class matrix<Value>;
  friend class const_matrix<Value>;
  friend class const_matrix_reference<Value>;
public:
  typedef matrix_reference<Value> self_type;
  typedef matrix_base<self_type> base_type;

  typedef typename base_type::size_type size_type;
  typedef typename base_type::value_type value_type;
  typedef typename base_type::reference reference;
  typedef typename base_type::const_reference const_reference;
  typedef typename base_type::pointer pointer;
  typedef typename base_type::const_pointer const_pointer;

  typedef typename base_type::vector vector;
  typedef typename base_type::const_vector const_vector;
  typedef typename base_type::vector_reference vector_reference;
  typedef typename base_type::const_vector_reference const_vector_reference;

  typedef typename base_type::matrix matrix;
  typedef typename base_type::const_matrix const_matrix;
  //typedef typename base_type::matrix_reference matrix_reference;
  typedef typename base_type::const_matrix_reference const_matrix_reference;

  typedef typename base_type::subscript subscript;
public:
  matrix_reference(pointer ptr, const size_type nrow,
		   const size_type ncol, const size_type ldim) :
    base_type(ptr, nrow, ncol, ldim)
  {
  }
  matrix_reference(const matrix& x)
  {
    base_type::alias(x);
  }
  matrix_reference(matrix&& x)
  {
    base_type::move(x);
  }
  matrix_reference(const_matrix&& x)
  {
    base_type::move(x);
  }
  matrix_reference(const matrix_reference& x)
  {
    base_type::alias(x);
  }
  matrix_reference(matrix_reference&& x)
  {
    base_type::move(x);
  }
  self_type& operator =(const matrix& x)
  {
    base_type::copy(x);

    return *this;
  }
  self_type& operator =(const const_matrix& x)
  {
    base_type::copy(x);

    return *this;
  }
  self_type& operator =(const matrix_reference& x)
  {
    base_type::copy(x);

    return *this;
  }
  self_type& operator =(const const_matrix_reference& x)
  {
    base_type::copy(x);

    return *this;
  }
  reference operator ()(const size_type row, const size_type col)
  {
    return base_type::at(row, col);
  }
  const_reference operator ()(const size_type row, const size_type col) const
  {
    return base_type::at(row, col);
  }
  matrix_reference operator ()(const subscript& row, const subscript& col)
  {
    return base_type::at(row, col);
  }
  const_matrix_reference operator ()(const subscript& row, const subscript& col) const
  {
    return base_type::at(row, col);
  }
  vector_reference row(const size_type pos)
  {
    base_type::row(pos);
  }
  const_vector_reference row(const size_type pos) const
  {
    base_type::row(pos);
  }
  vector_reference col(const size_type pos)
  {
    base_type::col(pos);
  }
  const_vector_reference col(const size_type pos) const
  {
    base_type::col(pos);
  }
  pointer ptr()
  {
    return base_type::ptr;
  }
  const_pointer ptr() const
  {
    return base_type::ptr_;
  }
  size_type nrow() const
  {
    return base_type::nrow_;
  }
  size_type ncol() const
  {
    return base_type::ncol_;
  }
  size_type ldim() const
  {
    return base_type::ldim_;
  }
};


template<typename Value>
class const_matrix_reference : private matrix_base<const_matrix_reference<Value> >
{
  friend class matrix<Value>;
  friend class const_matrix<Value>;
  friend class matrix_reference<Value>;
public:
  typedef const_matrix_reference<Value> self_type;
  typedef matrix_base<self_type> base_type;

  typedef typename base_type::size_type size_type;
  typedef typename base_type::value_type value_type;
  typedef typename base_type::reference reference;
  typedef typename base_type::const_reference const_reference;
  typedef typename base_type::pointer pointer;
  typedef typename base_type::const_pointer const_pointer;

  typedef typename base_type::vector vector;
  typedef typename base_type::const_vector const_vector;
  typedef typename base_type::vector_reference vector_reference;
  typedef typename base_type::const_vector_reference const_vector_reference;

  typedef typename base_type::matrix matrix;
  typedef typename base_type::const_matrix const_matrix;
  typedef typename base_type::matrix_reference matrix_reference;
  //typedef typename base_type::const_matrix_reference const_matrix_reference;

  typedef typename base_type::subscript subscript;
public:
  const_matrix_reference(pointer ptr, const size_type nrow,
			 const size_type ncol, const size_type ldim) :
    base_type(ptr, nrow, ncol, ldim)
  {
  }
  const_matrix_reference(const matrix& x)
  {
    base_type::alias(x);
  }
  const_matrix_reference(matrix&& x)
  {
    base_type::move(x);
  }
  const_matrix_reference(const const_matrix& x)
  {
    base_type::alias(x);
  }
  const_matrix_reference(const_matrix&& x)
  {
    base_type::move(x);
  }
  const_matrix_reference(const matrix_reference& x)
  {
    base_type::alias(x);
  }
  const_matrix_reference(matrix_reference&& x)
  {
    base_type::move(x);
  }
  const_matrix_reference(const const_matrix_reference& x)
  {
    base_type::alias(x);
  }
  const_matrix_reference(const_matrix_reference&& x)
  {
    base_type::move(x);
  }
  const_reference operator ()(const size_type row, const size_type col) const
  {
    return base_type::at(row, col);
  }
  const_matrix_reference operator ()(const subscript& row, const subscript& col) const
  {
    return base_type::at(row, col);
  }
  const_vector_reference row(const size_type pos) const
  {
    base_type::row(pos);
  }
  const_vector_reference col(const size_type pos) const
  {
    base_type::col(pos);
  }
  const_pointer ptr() const
  {
    return base_type::ptr;
  }
  size_type nrow() const
  {
    return base_type::nrow_;
  }
  size_type ncol() const
  {
    return base_type::ncol_;
  }
  size_type ldim() const
  {
    return base_type::ldim_;
  }
};


}  // namespace numeric
}  // namespace core
}  // namespace mcs


#endif
