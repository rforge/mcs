#ifndef _MCS_CORE_NUMERIC_MATRIX_HH_
#define _MCS_CORE_NUMERIC_MATRIX_HH_


#include <algorithm>
#include <cassert>
#include <utility>

#include "../../mcs.hh"

#include "subscript.hh"
#include "traits.hh"


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
  typedef typename traits<Derived>::reference_type reference_type;
  typedef typename traits<Derived>::const_reference_type const_reference_type;
  typedef typename traits<Derived>::pointer_type pointer_type;
  typedef typename traits<Derived>::const_pointer_type const_pointer_type;

  typedef typename traits<Derived>::vector_type vector_type;
  typedef typename traits<Derived>::const_vector_type const_vector_type;
  typedef typename traits<Derived>::vector_reference_type vector_reference_type;
  typedef typename traits<Derived>::const_vector_reference_type const_vector_reference_type;

  typedef typename traits<Derived>::matrix_type matrix_type;
  typedef typename traits<Derived>::const_matrix_type const_matrix_type;
  typedef typename traits<Derived>::matrix_reference_type matrix_reference_type;
  typedef typename traits<Derived>::const_matrix_reference_type const_matrix_reference_type;

  typedef typename traits<Derived>::subscript_type subscript_type;
protected:
  pointer_type buf_;
  pointer_type ptr_;
  size_type nrow_;
  size_type ncol_;
  size_type ldim_;
protected:
  vector_base() :
    buf_(nullptr),
    ptr_(nullptr),
    nrow_(0),
    ncol_(0),
    ldim_(0)
  {
  }
  explicit
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
  matrix_base(const pointer_type ptr, const size_type nrow,
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
  reference_type at(const size_type i, const size_type j)
  {
    MCS_ASSERT((0 <= i) && (i < nrow_), "invalid argument: i (matrix_base::at)");
    MCS_ASSERT((0 <= j) && (j < ncol_), "invalid argument: j (matrix_base::at)");

    return ptr_[j * ldim_ + i];
  }
  const_reference_type at(const size_type i, const size_type j) const
  {
    MCS_ASSERT((0 <= i) && (i < nrow_), "invalid argument: i (matrix_base::at)");
    MCS_ASSERT((0 <= j) && (j < ncol_), "invalid argument: j (matrix_base::at)");

    return ptr_[j * ldim_ + i];
  }
  matrix_reference_type at(const subscript_type& i, const subscript_type& j)
  {
    return matrix_reference_type(buf_ + j.pos * ldim_ + i.pos, i.len, j.len, ldim_);
  }
  const_matrix_reference_type at(const subscript_type& i, const subscript_type& j) const
  {
    return const_matrix_reference_type(buf_ + j.pos * ldim_ + i.pos, i.len, j.len, ldim_);
  }
  vector_reference_type row(const size_type i)
  {
    return vector_reference_type(buf_ + i, ncol_, ldim_);
  }
  const_vector_reference_type row(const size_type i) const
  {
    return const_vector_reference_type(buf_ + i, ncol_, ldim_);
  }
  vector_reference_type col(const size_type j)
  {
    return vector_reference_type(buf_ + j * ldim_, nrow_, 1);
  }
  const_vector_reference_type col(const size_type j) const
  {
    return const_vector_reference_type(buf_ + j * ldim_, nrow_, 1);
  }
  template<typename D>
  void copy(const matrix_base<D>& x)
  {
    MCS_ASSERT(nrow_ == x.nrow_, "invalid argument: x (matrix_base::copy)");
    MCS_ASSERT(ncol_ == x.ncol_, "invalid argument: x (matrix_base::copy)");

    pointer_type dst = buf_;
    const_pointer_type src = x.buf_;

    for (int j = 0; j < ncol_; ++j) {
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
  typedef typename base_type::value_type value_type;
  typedef typename base_type::reference_type reference_type;
  typedef typename base_type::const_reference_type const_reference_type;
  typedef typename base_type::pointer_type pointer_type;
  typedef typename base_type::const_pointer_type const_pointer_type;

  typedef typename base_type::vector_type vector_type;
  typedef typename base_type::const_vector_type const_vector_type;
  typedef typename base_type::vector_reference_type vector_reference_type;
  typedef typename base_type::const_vector_reference_type const_vector_reference_type;

  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::const_matrix_type const_matrix_type;
  typedef typename base_type::matrix_reference_type matrix_reference_type;
  typedef typename base_type::const_matrix_reference_type const_matrix_reference_type;

  typedef typename base_type::subscript_type subscript_type;
public:
  matrix(const size_type nrow, const size_type ncol) :
    base_type(nrow, ncol)
  {
  }
  matrix(const matrix_type& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  matrix(matrix_type&& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::move(x);
  }
  matrix(const const_matrix_type& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  matrix(const_matrix_type&& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::move(x);
  }
  matrix(const matrix_reference_type& x) :
    base_type(x.nrow, x.ncol_)
  {
    base_type::copy(x);
  }
  matrix(const const_matrix_reference_type& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  self_type& operator =(matrix_type x)
  {
    MCS_ASSERT(base_type::nrow_ == x.nrow_, "invalid argument: x (matrix::operator=)");
    MCS_ASSERT(base_type::ncol_ == x.ncol_, "invalid argument: x (matrix::operator=)");

    base_type::move(x);

    return *this;
  }
  reference_type operator ()(const size_type i, const size_type j)
  {
    return base_type::at(i, j);
  }
  const_reference_type operator ()(const size_type i, const size_type j) const
  {
    return base_type::at(i, j);
  }
  matrix_reference_type operator ()(const subscript_type& i, const subscript_type& j)
  {
    return base_type::at(i, j);
  }
  matrix_vector_reference_type operator ()(const subscript_type& i, const subscript_type& j) const
  {
    return base_type::at(i, j);
  }
  vector_reference_type row(const size_type i)
  {
    base_type::row(i);
  }
  const_vector_reference_type row(const size_type i) const
  {
    base_type::row(i);
  }
  vector_reference_type col(const size_type j)
  {
    base_type::col(j);
  }
  const_vector_reference_type col(const size_type j) const
  {
    base_type::col(j);
  }
  pointer_type ptr()
  {
    return base_type::ptr_;
  }
  const_pointer_type ptr() const
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
  typedef typename base_type::reference_type reference_type;
  typedef typename base_type::const_reference_type const_reference_type;
  typedef typename base_type::pointer_type pointer_type;
  typedef typename base_type::const_pointer_type const_pointer_type;

  typedef typename base_type::vector_type vector_type;
  typedef typename base_type::const_vector_type const_vector_type;
  typedef typename base_type::vector_reference_type vector_reference_type;
  typedef typename base_type::const_vector_reference_type const_vector_reference_type;

  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::const_matrix_type const_matrix_type;
  typedef typename base_type::matrix_reference_type matrix_reference_type;
  typedef typename base_type::const_matrix_reference_type const_matrix_reference_type;

  typedef typename base_type::subscript_type subscript_type;
public:
  explicit
  const_matrix(const size_type nrow, const size_type ncol) :
    base_type(nrow, ncol)
  {
  }
  const_matrix(const vector_type& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  const_matrix(matrix_type&& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::move(x);
  }
  const_matrix(const const_matrix_type& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  const_matrix(const_matrix_type&& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::move(x);
  }
  const_matrix(const matrix_reference_type& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  const_matrix(const const_matrix_reference_type& x) :
    base_type(x.nrow_, x.ncol_)
  {
    base_type::copy(x);
  }
  const_reference_type operator ()(const size_type i, const size_type j) const
  {
    return base_type::at(i, j);
  }
  const_matrix_reference_type operator ()(const subscript_type& i, const subscript_type& j) const
  {
    return base_type::at(i, j);
  }
  const_vector_reference_type row(const size_type i) const
  {
    base_type::row(i);
  }
  const_vector_reference_type col(const size_type j) const
  {
    base_type::col(j);
  }
  const_pointer_type ptr() const
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
  typedef typename base_type::reference_type reference_type;
  typedef typename base_type::const_reference_type const_reference_type;
  typedef typename base_type::pointer_type pointer_type;
  typedef typename base_type::const_pointer_type const_pointer_type;

  typedef typename base_type::vector_type vector_type;
  typedef typename base_type::const_vector_type const_vector_type;
  typedef typename base_type::vector_reference_type vector_reference_type;
  typedef typename base_type::const_vector_reference_type const_vector_reference_type;

  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::const_matrix_type const_matrix_type;
  typedef typename base_type::matrix_reference_type matrix_reference_type;
  typedef typename base_type::const_matrix_reference_type const_matrix_reference_type;

  typedef typename base_type::subscript_type subscript_type;
public:
  matrix_reference(const pointer_type ptr, const size_type nrow, const size_type ncol, const size_type ldim) :
    base_type(ptr, nrow, ncol, ldim)
  {
  }
  matrix_reference(const matrix_type& x)
  {
    base_type::alias(x);
  }
  matrix_reference(matrix_type&& x)
  {
    base_type::move(x);
  }
  matrix_reference(const_matrix_type&& x)
  {
    base_type::move(x);
  }
  matrix_reference(const matrix_reference_type& x)
  {
    base_type::alias(x);
  }
  matrix_reference(matrix_reference_type&& x)
  {
    base_type::move(x);
  }
  self_type& operator =(const matrix_type& x)
  {
    base_type::copy(x);

    return *this;
  }
  self_type& operator =(const const_matrix_type& x)
  {
    base_type::copy(x);

    return *this;
  }
  self_type& operator =(const matrix_reference_type& x)
  {
    base_type::copy(x);

    return *this;
  }
  self_type& operator =(const const_matrix_reference_type& x)
  {
    base_type::copy(x);

    return *this;
  }
  reference_type operator ()(size_type pos)
  {
    return base_type::at(pos);
  }
  const_reference_type operator ()(size_type pos) const
  {
    return base_type::at(pos);
  }
  matrix_reference_type operator ()(const subscript_type& i, const subscript_type& j)
  {
    return base_type::at(i, j);
  }
  const_matrix_reference_type operator ()(const subscript_type& i, const subscript_type& j) const
  {
    return base_type::at(i, j);
  }
  vector_reference_type row(const size_type i)
  {
    base_type::row(i);
  }
  const_vector_reference_type row(const size_type i) const
  {
    base_type::row(i);
  }
  vector_reference_type col(const size_type j)
  {
    base_type::col(j);
  }
  const_vector_reference_type col(const size_type j) const
  {
    base_type::col(j);
  }
  pointer_type ptr()
  {
    return base_type::ptr;
  }
  const_pointer_type ptr() const
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
  typedef typename base_type::reference_type reference_type;
  typedef typename base_type::const_reference_type const_reference_type;
  typedef typename base_type::pointer_type pointer_type;
  typedef typename base_type::const_pointer_type const_pointer_type;

  typedef typename base_type::vector_type vector_type;
  typedef typename base_type::const_vector_type const_vector_type;
  typedef typename base_type::vector_reference_type vector_reference_type;
  typedef typename base_type::const_vector_reference_type const_vector_reference_type;

  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::const_matrix_type const_matrix_type;
  typedef typename base_type::matrix_reference_type matrix_reference_type;
  typedef typename base_type::const_matrix_reference_type const_matrix_reference_type;

  typedef typename base_type::subscript_type subscript_type;
public:
  const_matrix_reference(const pointer_type ptr, const size_type len, const size_type inc) :
    base_type(ptr, len, inc)
  {
  }
  const_matrix_reference(const matrix_type& x)
  {
    base_type::alias(x);
  }
  const_matrix_reference(matrix_type&& x)
  {
    base_type::move(x);
  }
  const_matrix_reference(const const_matrix_type& x)
  {
    base_type::alias(x);
  }
  const_matrix_reference(const_matrix_type&& x)
  {
    base_type::move(x);
  }
  const_matrix_reference(const matrix_reference_type& x)
  {
    base_type::alias(x);
  }
  const_matrix_reference(matrix_reference_type&& x)
  {
    base_type::move(x);
  }
  const_matrix_reference(const const_matrix_reference_type& x)
  {
    base_type::alias(x);
  }
  const_matrix_reference(const_matrix_reference_type&& x)
  {
    base_type::move(x);
  }
  const_reference_type operator ()(const size_type i, const size_type j) const
  {
    return base_type::at(i, j);
  }
  const_matrix_reference_type operator ()(const subscript_type& i, const subscript_type& j) const
  {
    return base_type::at(i, j);
  }
  const_vector_reference_type row(const size_type i) const
  {
    base_type::row(i);
  }
  const_vector_reference_type col(const size_type j) const
  {
    base_type::col(j);
  }
  const_pointer_type ptr() const
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
