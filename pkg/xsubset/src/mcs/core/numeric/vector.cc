/**
 * @file vector.cc
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_NUMERIC_VECTOR_CC
#define MCS_CORE_NUMERIC_VECTOR_CC


#include "../../mcs.hh"

#include "range.hh"
#include "vector.hh"


#define VECTOR vector<Value, Alloc>
#define SUBVECTOR subvector<Value, Alloc>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
               typename Alloc>
      VECTOR::vector() :
        alloc_(),
        stor_(0),
        base_(0),
        inc_(0),
        len_(0)
      {
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(const size_type len) :
        alloc_(),
        stor_(alloc_.allocate(len)),
        base_(stor_),
        inc_(1),
        len_(len)
      {
        MCS_ASSERT(len >= 0);
      }


      template<typename Value,
               typename Alloc>
      template<typename V>
      VECTOR::vector(const size_type len,
                     const V val) :
        alloc_(),
        stor_(alloc_.allocate(len)),
        base_(stor_),
        inc_(1),
        len_(len)
      {
        MCS_ASSERT(len >= 0);

        fill(val);
      }


      template<typename Value,
               typename Alloc>
      template<typename V>
      VECTOR::vector(const size_type len,
                     const V* const data) :
        alloc_(),
        stor_(alloc_.allocate(len)),
        base_(stor_),
        inc_(1),
        len_(len)
      {
        MCS_ASSERT(len >= 0);

        fill(data);
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(const vector& vec) :
        alloc_(),
        stor_(alloc_.allocate(vec.len_)),
        base_(stor_),
        inc_(1),
        len_(vec.len_)
      {
        copy(vec);
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(vector&& vec) :
        alloc_(),
        stor_(0),
        base_(0),
        inc_(1),
        len_(0)
      {
        swap(vec);
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(subvector_type&& vec) :
        alloc_(),
        stor_(alloc_.allocate(vec.len_)),
        base_(stor_),
        inc_(1),
        len_(vec.len_)
      {
        copy(vec);
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(const pointer stor,
                     const size_type len) :
        alloc_(),
        stor_(stor),
        base_(stor_),
        inc_(1),
        len_(len)
      {
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(const pointer base,
                     const size_type inc,
                     const size_type len) :
        alloc_(),
        stor_(0),
        base_(base),
        inc_(inc),
        len_(len)
      {
      }


      template<typename Value,
               typename Alloc>
      VECTOR::~vector()
      {
        if (stor_)
          {
            alloc_.deallocate(stor_, len_);
          }
      }


      template<typename Value,
               typename Alloc>
      VECTOR&
      VECTOR::operator =(const vector& vec)
      {
        copy(vec);

        return *this;
      }


      template<typename Value,
               typename Alloc>
      VECTOR&
      VECTOR::operator =(vector&& vec)
      {
        swap(vec);

        return *this;
      }


      template<typename Value,
               typename Alloc>
      VECTOR&
      VECTOR::operator =(subvector_type&& vec)
      {
        copy(vec);

        return *this;
      }


      template<typename Value,
               typename Alloc>
      typename VECTOR::reference
      VECTOR::operator ()(const size_type pos)
      {
        MCS_ASSERT(pos <= len_);

        return at(pos);
      }


      template<typename Value,
               typename Alloc>
      typename VECTOR::const_reference
      VECTOR::operator ()(const size_type pos) const
      {
        MCS_ASSERT(pos <= len_);

        return at(pos);
      }


      template<typename Value,
               typename Alloc>
      typename VECTOR::subvector_type
      VECTOR::operator ()(const range_type rng)
      {
        const size_type pos = rng.pos();
        const size_type len = rng.open()? len_ - pos: rng.len();

        MCS_ASSERT((pos + len) <= len_);

        return subvector_type(&at(pos), inc_, len);
      }


      template<typename Value,
               typename Alloc>
      const typename VECTOR::subvector_type
      VECTOR::operator ()(const range_type rng) const
      {
        const size_type pos = rng.pos();
        const size_type len = rng.open()? len_ - pos: rng.len();

        MCS_ASSERT((pos + len) <= len_);

        return subvector_type(&at(pos), inc_, len);
      }


      template<typename Value,
               typename Alloc>
      typename VECTOR::size_type
      VECTOR::len() const
      {
        return len_;
      }


      template<typename Value,
               typename Alloc>
      typename VECTOR::size_type
      VECTOR::inc() const
      {
        return inc_;
      }


      template<typename Value,
               typename Alloc>
      template<typename V>
      void
      VECTOR::fill(const V val)
      {
	pointer ptr = base_;
        for (size_type i = 0; i < len_; ++i)
          {
            *ptr = static_cast<value_type>(val);
	    ptr += inc_;
          }
      }


      template<typename Value,
	       typename Alloc>
      template<typename V>
      void
      VECTOR::fill(const V* data)
      {
	pointer ptr = base_;
	for (size_type i = 0; i < len_; ++i)
	  {
            ptr = static_cast<value_type>(*data++);
            ptr += inc_;
	  }
      }


      template<typename Value,
               typename Alloc>
      void
      VECTOR::copy(const vector& vec)
      {
        MCS_ASSERT(len_ == vec.len_);

	const_pointer src = vec.base_;
	pointer dst = base_;
        for (size_type i = 0; i < len_; ++i)
          {
	    *dst = *src;
	    src += vec.inc_;
	    dst += inc_;
          }
      }


      template<typename Value,
               typename Alloc>
      void
      VECTOR::swap(vector& vec)
      {
	// TODO: swap allocator
        std::swap(stor_, vec.stor_);
        std::swap(base_, vec.base_);
        std::swap(inc_, vec.inc_);
        std::swap(len_, vec.len_);
      }


      template<typename Value,
               typename Alloc>
      typename VECTOR::reference
      VECTOR::at(const size_type pos) const
      {
        return base_[pos * inc_];
      }


      template<typename Value,
               typename Alloc>
      SUBVECTOR::subvector(pointer base,
                           size_type inc,
                           size_type len)
        : vector_type(base, inc, len)
      {
      }


      template<typename Value,
               typename Alloc>
      SUBVECTOR::~subvector()
      {
      }


      template<typename Value,
               typename Alloc>
      SUBVECTOR&
      SUBVECTOR::operator =(const vector_type& vec)
      {
        MCS_ASSERT(len_ == vec.len_);

        copy(vec);

        return *this;
      }


      template<typename Value,
               typename Alloc>
      SUBVECTOR&
      SUBVECTOR::operator =(vector_type&& vec)
      {
        MCS_ASSERT(len_ == vec.len_);

        copy(vec);

        return *this;
      }


    }

  }

}


#undef VECTOR
#undef SUBVECTOR


#define VECTOR mcs::core::numeric::vector<Value, Alloc>
#define SUBVECTOR mcs::core::numeric::subvector<Value, Alloc>


namespace std
{


  template<typename Value,
	   typename Alloc>
  void
  swap(VECTOR& x, VECTOR& y)
  {
    x.swap(y);
  }


}


#undef VECTOR
#undef SUBVECTOR


#endif
