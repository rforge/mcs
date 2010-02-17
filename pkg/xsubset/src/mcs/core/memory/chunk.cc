
//          Copyright Marc Hofmann 2009 - 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file chunk.cc
 *
 * @author Marc Hofmann
 *
 * Work in progress.
 *
 * Derived from Boost Pool.
 */

#ifndef MCS_CORE_MEMORY_CHUNK_CC
#define MCS_CORE_MEMORY_CHUNK_CC


#include "chunk.hh"


#define CHUNK chunk<Size>


namespace mcs
{

  namespace core
  {

    namespace memory
    {


      template<typename Size>
      CHUNK::chunk(char* const ptr) : ptr_(ptr)
      {
      }


      template<typename Size>
      CHUNK::operator char*() const
      {
        return ptr_;
      }


      template<typename Size>
      CHUNK&
      CHUNK::operator =(char* const ptr)
      {
        ptr_ = ptr;
        return *this;
      }


      template<typename Size>
      CHUNK
      CHUNK::operator +(const CHUNK::size_type dist)
      {
        return chunk(ptr_ + dist);
      }


      template<typename Size>
      CHUNK
      CHUNK::operator -(const CHUNK::size_type dist)
      {
        return chunk(ptr_ - dist);
      }


      template<typename Size>
      CHUNK&
      CHUNK::operator +=(const CHUNK::size_type dist)
      {
        ptr_ += dist;
        return *this;
      }


      template<typename Size>
      CHUNK&
      CHUNK::operator -=(const CHUNK::size_type dist)
      {
        ptr_ -= dist;
        return *this;
      }


      template<typename Size>
      bool
      CHUNK::operator ==(char* const ptr)
      {
        return ptr_ == ptr;
      }


      template<typename Size>
      bool
      CHUNK::operator !=(char* const ptr)
      {
        return ptr_ != ptr;
      }


      template<typename Size>
      bool
      CHUNK::operator >(char* const ptr)
      {
        return ptr_ > ptr;
      }


      template<typename Size>
      CHUNK
      CHUNK::next() const
      {
        return link::next(*this);
      }


      template<typename Size>
      void
      CHUNK::next(const CHUNK next) const
      {
        link::next(*this) = next.ptr_;
      }


      template<typename Size>
      char*&
      CHUNK::link::next(const CHUNK chunk)
      {
        return reinterpret_cast<link*>(chunk.ptr_)->next_;
      }


    }

  }

}


#endif
