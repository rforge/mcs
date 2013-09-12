#ifndef _MCS_CORE_NUMERIC_VECTOR_HH_
#define _MCS_CORE_NUMERIC_VECTOR_HH_


#include <cassert>
#include <utility>

#include "../../mcs.hh"

#include "subscript.hh"
#include "traits.hh"


namespace mcs     {
namespace core    {
namespace numeric {


// vector base
template<typename Derived>
class vector_base
{
  template<typename D>
  friend class vector_base;
public:
  typedef vector_base<Derived> self_type;
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

  typedef typename traits<Derived>::subscript_type subscript_type;
protected:
  pointer_type buf_;
  pointer_type ptr_;
  size_type len_;
  size_type inc_;
protected:
  vector_base() :
    buf_(nullptr),
    ptr_(nullptr),
    len_(0),
    inc_(0)
  {
  }
  explicit
  vector_base(const size_type len) :
    buf_(new value_type[len]),
    ptr_(buf_),
    len_(len),
    inc_(1)
  {
    MCS_ASSERT(len >= 0, "invalid argument: len (vector_base::vector_base)");
  }
  vector_base(const pointer_type ptr, const size_type len, const size_type inc) :
    buf_(nullptr),
    ptr_(ptr),
    len_(len),
    inc_(inc)
  {
    MCS_ASSERT(ptr != nullptr, "invalid argument: ptr (vector_base::vector_base)");
    MCS_ASSERT(len >= 0, "invalid argument: len (vector_base::vector_base)");
    MCS_ASSERT(inc >= 0, "invalid argument: inc (vector_base::vector_base)");
  }
  ~vector_base()
  {
    if (buf_ != nullptr)
      {
        delete [] buf_;
        buf_ = nullptr;
      }
  }
  reference_type at(const size_type pos)
  {
    MCS_ASSERT((0 <= pos) && (pos < len_), "invalid argument: pos (vector_base::at)");

    return ptr_[pos * inc_];
  }
  const_reference_type at(const size_type pos) const
  {
    MCS_ASSERT((0 <= pos) && (pos < len_), "invalid argument: pos (vector_base::at)");

    return ptr_[pos * inc_];
  }
  vector_reference_type at(const subscript_type& pos)
  {
    return vector_reference_type(ptr_ + pos.pos * inc_, pos.len, pos.inc * inc_);
  }
  const_vector_reference_type at(const subscript_type& pos) const
  {
    return const_vector_reference_type(ptr_ + pos.pos * inc_, pos.len, pos.inc * inc_);
  }
  template<typename D>
  void copy(const vector_base<D>& x)
  {
    MCS_ASSERT(len_ == x.len_, "invalid argument: x (vector_base::copy)");

    const pointer_type end = ptr_ + len_ * inc_;

    pointer_type dst = ptr_;
    const_pointer_type src = x.ptr_;

    while (dst < end)
      {
        *dst = *src;

        dst += inc_;
        src += x.inc_;
      }
  }
  template<typename D>
  void move(vector_base<D>& x)
  {
    std::swap(buf_, x.buf_);
    std::swap(ptr_, x.ptr_);
    std::swap(len_, x.len_);
    std::swap(inc_, x.inc_);
  }
  template<typename D>
  void alias(const vector_base<D>& x)
  {
    ptr_ = x.ptr_;
    len_ = x.len_;
    inc_ = x.inc_;
  }
};


// vector, const_vector, vector_reference, const_vector_reference
template<typename Value>
class vector : private vector_base<vector<Value> >
{
  friend class const_vector<Value>;
  friend class vector_reference<Value>;
  friend class const_vector_reference<Value>;
public:
  typedef vector<Value> self_type;
  typedef vector_base<self_type> base_type;

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

  typedef typename base_type::subscript_type subscript_type;
public:
  explicit
  vector(const size_type len) :
    base_type(len)
  {
  }
  vector(const vector_type& x) :
    base_type(x.len_)
  {
    base_type::copy(x);
  }
  vector(vector_type&& x) :
    base_type(x.len_)
  {
    base_type::move(x);
  }
  vector(const const_vector_type& x) :
    base_type(x.len_)
  {
    base_type::copy(x);
  }
  vector(const_vector_type&& x) :
    base_type(x.len_)
  {
    base_type::move(x);
  }
  vector(const vector_reference_type& x) :
    base_type(x.len_)
  {
    base_type::copy(x);
  }
  vector(const const_vector_reference_type& x) :
    base_type(x.len_)
  {
    base_type::copy(x);
  }
  self_type& operator =(vector_type x)
  {
    MCS_ASSERT(base_type::len_ == x.len_, "invalid argument: x (vector::operator=)");

    base_type::move(x);

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
  vector_reference_type operator ()(const subscript_type& pos)
  {
    return base_type::at(pos);
  }
  const_vector_reference_type operator ()(const subscript_type& pos) const
  {
    return base_type::at(pos);
  }
  pointer_type ptr()
  {
    return base_type::ptr_;
  }
  const_pointer_type ptr() const
  {
    return base_type::ptr_;
  }
  size_type len() const
  {
    return base_type::len_;
  }
  size_type inc() const
  {
    return base_type::inc_;
  }
};


template<typename Value>
class const_vector : private vector_base<const_vector<Value> >
{
  friend class vector<Value>;
  friend class vector_reference<Value>;
  friend class const_vector_reference<Value>;
public:
  typedef const_vector<Value> self_type;
  typedef vector_base<self_type> base_type;

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

  typedef typename base_type::subscript_type subscript_type;
public:
  explicit
  const_vector(const size_type len) :
    base_type(len)
  {
  }
  const_vector(const vector_type& x) :
    base_type(x.len_)
  {
    base_type::copy(x);
  }
  const_vector(vector_type&& x) :
    base_type(x.len_)
  {
    base_type::move(x);
  }
  const_vector(const const_vector_type& x) :
    base_type(x.len_)
  {
    base_type::copy(x);
  }
  const_vector(const_vector_type&& x) :
    base_type(x.len_)
  {
    base_type::move(x);
  }
  const_vector(const vector_reference_type& x) :
    base_type(x.len_)
  {
    base_type::copy(x);
  }
  const_vector(const const_vector_reference_type& x) :
    base_type(x.len_)
  {
    base_type::copy(x);
  }
  const_reference_type operator ()(size_type pos) const
  {
    return base_type::at(pos);
  }
  const_vector_reference_type operator ()(const subscript_type& pos) const
  {
    return base_type::at(pos);
  }
  const_pointer_type ptr() const
  {
    return base_type::ptr_;
  }
  size_type len() const
  {
    return base_type::len_;
  }
  size_type inc() const
  {
    return base_type::inc_;
  }
};


template<typename Value>
class vector_reference : private vector_base<vector_reference<Value> >
{
  friend class vector<Value>;
  friend class const_vector<Value>;
  friend class const_vector_reference<Value>;
public:
  typedef vector_reference<Value> self_type;
  typedef vector_base<self_type> base_type;

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

  typedef typename base_type::subscript_type subscript_type;
public:
  vector_reference(const pointer_type ptr, const size_type len, const size_type inc) :
    base_type(ptr, len, inc)
  {
  }
  vector_reference(const vector_type& x)
  {
    base_type::alias(x);
  }
  vector_reference(vector_type&& x)
  {
    base_type::move(x);
  }
  vector_reference(const_vector_type&& x)
  {
    base_type::move(x);
  }
  vector_reference(const vector_reference_type& x)
  {
    base_type::alias(x);
  }
  vector_reference(vector_reference_type&& x)
  {
    base_type::move(x);
  }
  self_type& operator =(const vector_type& x)
  {
    base_type::copy(x);

    return *this;
  }
  self_type& operator =(const const_vector_type& x)
  {
    base_type::copy(x);

    return *this;
  }
  self_type& operator =(const vector_reference_type& x)
  {
    base_type::copy(x);

    return *this;
  }
  self_type& operator =(const const_vector_reference_type& x)
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
  vector_reference_type operator ()(const subscript_type& pos)
  {
    return base_type::at(pos);
  }
  const_vector_reference_type operator ()(const subscript_type& pos) const
  {
    return base_type::at(pos);
  }
  pointer_type ptr()
  {
    return base_type::ptr_;
  }
  const_pointer_type ptr() const
  {
    return base_type::ptr_;
  }
  size_type len() const
  {
    return base_type::len_;
  }
  size_type inc() const
  {
    return base_type::inc_;
  }
};


template<typename Value>
class const_vector_reference : private vector_base<const_vector_reference<Value> >
{
  friend class vector<Value>;
  friend class const_vector<Value>;
  friend class vector_reference<Value>;
public:
  typedef const_vector_reference<Value> self_type;
  typedef vector_base<self_type> base_type;

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

  typedef typename base_type::subscript_type subscript_type;
public:
  const_vector_reference(const pointer_type ptr, const size_type len, const size_type inc) :
    base_type(ptr, len, inc)
  {
  }
  const_vector_reference(const vector_type& x)
  {
    base_type::alias(x);
  }
  const_vector_reference(vector_type&& x)
  {
    base_type::move(x);
  }
  const_vector_reference(const const_vector_type& x)
  {
    base_type::alias(x);
  }
  const_vector_reference(const_vector_type&& x)
  {
    base_type::move(x);
  }
  const_vector_reference(const vector_reference_type& x)
  {
    base_type::alias(x);
  }
  const_vector_reference(vector_reference_type&& x)
  {
    base_type::move(x);
  }
  const_vector_reference(const const_vector_reference_type& x)
  {
    base_type::alias(x);
  }
  const_vector_reference(const_vector_reference_type&& x)
  {
    base_type::move(x);
  }
  const_reference_type operator ()(size_type pos) const
  {
    return base_type::at(pos);
  }
  const_vector_reference_type operator ()(const subscript_type& pos) const
  {
    return base_type::at(pos);
  }
  const_pointer_type ptr() const
  {
    return base_type::ptr_;
  }
  size_type len() const
  {
    return base_type::len_;
  }
  size_type inc() const
  {
    return base_type::inc_;
  }
};


}  // namespace numeric
}  // namespace core
}  // namespace mcs


#endif
