
//          Copyright Marc Hofmann 2009 - 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file chunk.hh
 *
 * @author Marc Hofmann
 *
 * Work in progress.
 *
 * Derived from Boost Pool.
 */

#ifndef MCS_CORE_MEMORY_CHUNK_HH
#define MCS_CORE_MEMORY_CHUNK_HH


namespace mcs
{

  namespace core
  {

    namespace memory
    {


      template<typename Size>
      class chunk
      {

      public:

        typedef Size size_type;

        typedef chunk<size_type> self_type;


      private:

        char* ptr_;


      public:

        chunk(char* ptr = 0);

        operator char*() const;

        self_type&
        operator =(char* ptr);

        self_type
        operator +(size_type diff);

        self_type
        operator -(size_type diff);

        self_type&
        operator +=(size_type diff);

        self_type&
        operator -=(size_type diff);

        bool
        operator ==(char* ptr);

        bool
        operator !=(char* ptr);

        bool
        operator >(char* ptr);

        self_type
        next() const;

        void
        next(self_type next) const;


      private:

        struct link
        {
          char* next_;

          static char*& next(self_type chunk);
        };

      };


    }

  }

}


#include "chunk.cc"
#endif
