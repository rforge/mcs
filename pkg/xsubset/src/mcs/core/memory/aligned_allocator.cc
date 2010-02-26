/**
 * @file aligned_allocator.cc
 *
 * @author Marc Hofmann
 *
 * Work in progress.
 */

#ifndef MCS_CORE_MEMORY_ALIGNED_ALLOCATOR_CC
#define MCS_CORE_MEMORY_ALIGNED_ALLOCATOR_CC


#define ALIGNED_ALLOC aligned_allocator<Value, Align>


namespace mcs
{

  namespace core
  {

    namespace memory
    {


      template<typename Value,
	       size_t Align>
      ALIGNED_ALLOC::aligned_allocator()
      {
      }


      template<typename Value,
	       size_t Align>
      ALIGNED_ALLOC::aligned_allocator(const aligned_allocator& alloc)
      {
      }


      template<typename Value,
	       size_t Align>
      template<typename U>
      ALIGNED_ALLOC::aligned_allocator(const aligned_allocator<U, Align>& alloc)
      {
      }


      template<typename Value,
	       size_t Align>
      ALIGNED_ALLOC::~aligned_allocator()
      {
      }


      template<typename Value,
	       size_t Align>
      template<typename U>
      ALIGNED_ALLOC&
      ALIGNED_ALLOC::operator =(const aligned_allocator<U, Align>& alloc)
      {
	return *this;
      }


      template<typename Value,
	       size_t Align>
      typename ALIGNED_ALLOC::pointer
      ALIGNED_ALLOC::address(reference ref) const
      {
	return &ref;
      }


      template<typename Value,
	       size_t Align>
      typename ALIGNED_ALLOC::const_pointer
      ALIGNED_ALLOC::address(const_reference ref) const
      {
	return &ref;
      }

      template<typename Value,
	       size_t Align>
      typename ALIGNED_ALLOC::pointer
      ALIGNED_ALLOC::allocate(const size_type num,
			      const void*)
      {
	void* mem = malloc(num * sizeof(value_type) + ALIGNMENT + sizeof(void*));
	void* ptr = (void*) (((uintptr_t) mem + sizeof(void*) + ALIGNMENT) & MASK);
	((void**) ptr)[-1] = mem;
	return (pointer) ptr;
      }


      template<typename Value,
	       size_t Align>
      void
      ALIGNED_ALLOC::deallocate(const pointer ptr,
				const size_type num)
      {
	free(((void**) ptr)[-1]);
      }


      template<typename Value,
	       size_t Align>
      void
      ALIGNED_ALLOC::construct(const pointer ptr,
			       const_reference val)
      {
	new ((void*) ptr) value_type(val);
      }


      template<typename Value,
	       size_t Align>
      void
      ALIGNED_ALLOC::destroy(const pointer ptr)
      {
	ptr->~value_type();
      }


    }

  }

}


#undef ALIGNED_ALLOC


#endif
