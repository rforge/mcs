
//          Copyright Marc Hofmann 2009 - 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file block.cc
 *
 * @author Marc Hofmann
 *
 * Work in progress.
 *
 * Derived from Boost Pool.
 */

#ifndef MCS_CORE_MEMORY_BLOCK_CC
#define MCS_CORE_MEMORY_BLOCK_CC


#include "block.hh"


#define BLOCK block<Size>


namespace mcs
{

  namespace core
  {

    namespace memory
    {


      template<typename Size>
      BLOCK::block(char* const ptr) : ptr_(ptr)
      {
      }


      template<typename Size>
      BLOCK::block(char* const ptr,
                   const typename BLOCK::size_type size) :
        ptr_(ptr)
      {
        descriptor::size(*this) = size;
      }


      template<typename Size>
      char*
      BLOCK::ptr() const
      {
        return ptr_;
      }


      template<typename Size>
      typename BLOCK::size_type
      BLOCK::size() const
      {
        return descriptor::size(*this);
      }


      template<typename Size>
      typename BLOCK::size_type
      BLOCK::chunk_area() const
      {
        return size() - sizeof(descriptor);
      }


      template<typename Size>
      BLOCK
      BLOCK::next() const
      {
        return block(descriptor::next_ptr(*this));
      }


      template<typename Size>
      void
      BLOCK::next(const BLOCK next)
      {
        descriptor::next_ptr(*this) = next.ptr_;
      }


      template<typename Size>
      typename BLOCK::chunk_type
      BLOCK::begin() const
      {
        return chunk_type(ptr_ + sizeof(descriptor));
      }


      template<typename Size>
      typename BLOCK::chunk_type
      BLOCK::end() const
      {
        return chunk_type(ptr_ + sizeof(descriptor) + chunk_area());
      }


      template<typename Size>
      typename BLOCK::size_type&
      BLOCK::descriptor::size(const BLOCK block)
      {
        return reinterpret_cast<descriptor*>(block.ptr_)->size_;
      }


      template<typename Size>
      char*&
      BLOCK::descriptor::next_ptr(const BLOCK block)
      {
        return reinterpret_cast<descriptor*>(block.ptr_)->next_ptr_;
      }


    }

  }

}


#endif
