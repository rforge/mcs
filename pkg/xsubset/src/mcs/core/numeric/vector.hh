#ifndef _MCS_CORE_NUMERIC_VECTOR_HH_
#define _MCS_CORE_NUMERIC_VECTOR_HH_


#include <cassert>
#include <cstddef>
#include <utility>

#include "subscript.hh"


namespace mcs     {
namespace core    {
namespace numeric {


// forward declarations
template<typename Value>
class vector;

template<typename Value>
class const_vector;

template<typename Value>
class vector_reference;

template<typename Value>
class const_vector_reference;


// vector traits
template<typename Value>
struct vector_traits
{
  typedef std::size_t size_type;
  typedef Value value_type;
  typedef Value& reference_type;
  typedef const Value& const_reference_type;
  typedef Value* pointer_type;
  typedef const Value* const_pointer_type;

  typedef vector<Value> vector_type;
  typedef const_vector<Value> const_vector_type;
  typedef vector_reference<Value> vector_reference_type;
  typedef const_vector_reference<Value> const_vector_reference_type;

  typedef subscript<size_type> subscript_type;
};

template<typename Value>
struct vector_traits<vector<Value> > : vector_traits<Value>
{
};

template<typename Value>
struct vector_traits<const_vector<Value> > : vector_traits<Value>
{
};

template<typename Value>
struct vector_traits<vector_reference<Value> > : vector_traits<Value>
{
};

template<typename Value>
struct vector_traits<const_vector_reference<Value> > : vector_traits<Value>
{
};



// vector base
template<typename Derived>
class vector_base
{
  template<typename D>
  friend class vector_base;
public:
  typedef vector_base<Derived> self_type;
  typedef Derived derived_type;

  typedef typename vector_traits<Derived>::size_type size_type;
  typedef typename vector_traits<Derived>::value_type value_type;
  typedef typename vector_traits<Derived>::reference_type reference_type;
  typedef typename vector_traits<Derived>::const_reference_type const_reference_type;
  typedef typename vector_traits<Derived>::pointer_type pointer_type;
  typedef typename vector_traits<Derived>::const_pointer_type const_pointer_type;

  typedef typename vector_traits<Derived>::vector_type vector_type;
  typedef typename vector_traits<Derived>::const_vector_type const_vector_type;
  typedef typename vector_traits<Derived>::vector_reference_type vector_reference_type;
  typedef typename vector_traits<Derived>::const_vector_reference_type const_vector_reference_type;

  typedef typename vector_traits<Derived>::subscript_type subscript_type;
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
  }
  vector_base(const pointer_type ptr, const size_type len, const size_type inc) :
    buf_(nullptr),
    ptr_(ptr),
    len_(len),
    inc_(inc)
  {
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
    assert((0 <= pos) && (pos < len_));

    return ptr_[pos * inc_];
  }
  const_reference_type at(const size_type pos) const
  {
    assert((0 <= pos) && (pos < len_));

    return ptr_[pos * inc_];
  }
  vector_reference_type at(const subscript_type& slice)
  {
    return vector_reference_type(ptr_ + slice.pos * inc_, slice.len, slice.inc * inc_);
  }
  const_vector_reference_type at(const subscript_type& slice) const
  {
    return const_vector_reference_type(ptr_ + slice.pos * inc_, slice.len, slice.inc * inc_);
  }
  template<typename D>
  void copy(const vector_base<D>& x)
  {
    assert(len_ == x.len_);

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
    assert(base_type::len_ == x.len_);

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
  vector_reference_type operator ()(const subscript_type& slice)
  {
    return base_type::at(slice);
  }
  const_vector_reference_type operator ()(const subscript_type& slice) const
  {
    return base_type::at(slice);
  }
  pointer_type ptr()
  {
    return base_type::ptr_ + base_type::off_;
  }
  const_pointer_type ptr() const
  {
    return base_type::ptr_ + base_type::off_;
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
  const_vector_reference_type operator ()(const subscript_type& slice) const
  {
    return base_type::at(slice);
  }
  const_pointer_type ptr() const
  {
    return base_type::ptr_ + base_type::off_;
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
  vector_reference_type operator ()(const subscript_type& slice)
  {
    return base_type::at(slice);
  }
  const_vector_reference_type operator ()(const subscript_type& slice) const
  {
    return base_type::at(slice);
  }
  pointer_type ptr()
  {
    return base_type::ptr_ + base_type::off_;
  }
  const_pointer_type ptr() const
  {
    return base_type::ptr_ + base_type::off_;
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
  const_vector_reference_type operator ()(const subscript_type& slice) const
  {
    return base_type::at(slice);
  }
  const_pointer_type ptr() const
  {
    return base_type::ptr_ + base_type::off_;
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
