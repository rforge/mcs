/**
 * @file aligned_allocator.hh
 *
 * @author Marc Hofmann
 *
 * Work in progress.
 */

#ifndef MCS_CORE_MEMORY_ALIGNED_ALLOCATOR_HH
#define MCS_CORE_MEMORY_ALIGNED_ALLOCATOR_HH


namespace mcs
{

  namespace core
  {

    namespace memory
    {


      template<typename Value,
	       size_t Align>
      class aligned_allocator
      {

      private:

	static const size_t ALIGNMENT = Align;

	static const uintptr_t MASK = ~(uintptr_t) (ALIGNMENT - 1);


      public:

	typedef Value value_type;

	typedef Value* pointer;

	typedef const Value* const_pointer;

	typedef Value& reference;

	typedef const Value& const_reference;

	typedef size_t size_type;

	typedef ptrdiff_t difference_type;


	template <class U>
	struct rebind {

	  typedef aligned_allocator<U, ALIGNMENT> other;

	};


      public:

	aligned_allocator();

	aligned_allocator(const aligned_allocator& alloc);

	template<typename U>
	aligned_allocator(const aligned_allocator<U, Align>& alloc);

	~aligned_allocator();


	template<typename U>
	aligned_allocator&
	operator =(const aligned_allocator<U, Align>& alloc);

	pointer
	address(reference ref) const;

	const_pointer
	address(const_reference ref) const;

	pointer
	allocate(size_type num,
		 const void* = 0);

	void
	deallocate(pointer ptr,
		   size_type num);

	void
	construct(pointer ptr,
		  const_reference val);

	void
	destroy(pointer ptr);

      };


    }

  }

}


#include "aligned_allocator.cc"
#endif
