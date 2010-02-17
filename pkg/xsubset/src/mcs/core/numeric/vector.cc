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
      VECTOR::vector_impl::vector_impl() :
        allocator_type(),
        stor(0),
        base(0),
        inc(0),
        len(0)
      {
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector_impl::vector_impl(const allocator_type& alloc) :
        allocator_type(alloc),
        stor(0),
        base(0),
        inc(0),
        len(0)
      {
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector() :
        impl_()
      {
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(const size_type len) :
        impl_()
      {
        MCS_ASSERT(len >= 0);

        impl_.stor = impl_.allocate(len);
        impl_.base = impl_.stor;
        impl_.inc = 1;
        impl_.len = len;
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(const typename VECTOR::size_type len,
                     const typename VECTOR::value_type val) :
        impl_()
      {
        MCS_ASSERT(len >= 0);

        impl_.stor = impl_.allocate(len);
        impl_.base = impl_.stor;
        impl_.inc = 1;
        impl_.len = len;

        fill(val);
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(const VECTOR& vec) :
        impl_()
      {
        if (vec.impl_.base)
          {
            impl_.stor = impl_.allocate(vec.impl_.len);
            impl_.base = impl_.stor;
            impl_.inc = 1;
            impl_.len = vec.impl_.len;

            copy(vec);
          }
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(VECTOR&& vec) :
        impl_()
      {
        swap(vec);
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(typename VECTOR::subvector_type&& vec) :
        impl_()
      {
        impl_.stor = impl_.allocate(vec.impl_.len);
        impl_.base = impl_.stor;
        impl_.inc = 1;
        impl_.len = vec.impl_.len;

        copy(vec);
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(const typename VECTOR::pointer stor,
                     const typename VECTOR::size_type len) :
        impl_()
      {
        impl_.stor = stor;
        impl_.base = impl_.stor;
        impl_.inc = 1;
        impl_.len = len;
      }


      template<typename Value,
               typename Alloc>
      VECTOR::vector(const typename VECTOR::pointer base,
                     const typename VECTOR::size_type inc,
                     const typename VECTOR::size_type len) :
        impl_()
      {
        impl_.stor = 0;
        impl_.base = base;
        impl_.inc = inc;
        impl_.len = len;
      }


      template<typename Value,
               typename Alloc>
      VECTOR::~vector()
      {
        if (impl_.stor)
          {
            impl_.deallocate(impl_.stor, 0);
          }
      }


      template<typename Value,
               typename Alloc>
      VECTOR&
      VECTOR::operator =(const VECTOR vec)
      {
        swap(vec);

        return *this;
      }


      template<typename Value,
               typename Alloc>
      typename VECTOR::reference
      VECTOR::operator ()(const typename VECTOR::size_type pos)
      {
        MCS_ASSERT(pos <= impl_.len);

        return at(pos);
      }


      template<typename Value,
               typename Alloc>
      typename VECTOR::const_reference
      VECTOR::operator ()(const typename VECTOR::size_type pos) const
      {
        MCS_ASSERT(pos <= impl_.len);

        return at(pos);
      }


      template<typename Value,
               typename Alloc>
      typename VECTOR::subvector_type
      VECTOR::operator ()(const typename VECTOR::range_type rng)
      {
        const size_type pos = rng.pos();
        const size_type len = rng.open()? impl_.len - pos: rng.len();

        MCS_ASSERT((pos + len) <= impl_.len);

        return subvector_type(&at(pos), impl_.inc, len);
      }


      template<typename Value,
               typename Alloc>
      const typename VECTOR::subvector_type
      VECTOR::operator ()(const typename VECTOR::range_type rng) const
      {
        const size_type pos = rng.pos();
        const size_type len = rng.open()? impl_.len - pos: rng.len();

        MCS_ASSERT((pos + len) <= impl_.len);

        return subvector_type(&at(pos), impl_.inc, len);
      }


      template<typename Value,
               typename Alloc>
      typename VECTOR::size_type
      VECTOR::len() const
      {
        return impl_.len;
      }


      template<typename Value,
               typename Alloc>
      typename VECTOR::size_type
      VECTOR::inc() const
      {
        return impl_.inc;
      }


      template<typename Value,
               typename Alloc>
      void
      VECTOR::fill(const typename VECTOR::value_type val)
      {
	const size_type len = impl_.len;

	pointer ptr = impl_.base;
        for (size_type i = 0; i < len; ++i)
          {
            *ptr = val;
	    ptr += impl_.inc;
          }
      }


      template<typename Value,
               typename Alloc>
      void
      VECTOR::copy(const VECTOR& vec)
      {
        MCS_ASSERT(impl_.len == vec.impl_.len);

	const size_type isrc = vec.impl_.inc;
	const size_type idst = impl_.inc;
	const size_type len = impl_.len;

	const_pointer src = vec.impl_.base;
	pointer dst = impl_.base;
        for (size_type i = 0; i < len; ++i)
          {
	    *dst = *src;
	    src += isrc;
	    dst += idst;
          }
      }


      template<typename Value,
               typename Alloc>
      void
      VECTOR::swap(VECTOR& vec)
      {
        // TODO: swap allocator
        std::swap(impl_.stor, vec.impl_.stor);
        std::swap(impl_.base, vec.impl_.base);
        std::swap(impl_.inc, vec.impl_.inc);
        std::swap(impl_.len, vec.impl_.len);
      }


      template<typename Value,
               typename Alloc>
      typename VECTOR::reference
      VECTOR::at(const typename VECTOR::size_type pos) const
      {
        return impl_.base[pos * impl_.inc];
      }


      template<typename Value,
               typename Alloc>
      SUBVECTOR::subvector(const typename SUBVECTOR::pointer base,
                           const typename SUBVECTOR::size_type inc,
                           const typename SUBVECTOR::size_type len)
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
      SUBVECTOR::operator =(const typename SUBVECTOR::vector_type& vec)
      {
        MCS_ASSERT(this->impl_.len == vec.impl_.len);

        vector_type::copy(vec);

        return *this;
      }


      template<typename Value,
               typename Alloc>
      SUBVECTOR&
      SUBVECTOR::operator =(typename SUBVECTOR::vector_type&& vec)
      {
        MCS_ASSERT(this->impl_.len == vec.impl_.len);

        vector_type::copy(vec);

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
