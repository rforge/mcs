
//          Copyright Marc Hofmann 2009 - 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file block.hh
 *
 * @author Marc Hofmann
 *
 * Work in progress.
 *
 * Derived from Boost Pool.
 */

#ifndef MCS_CORE_MEMORY_BLOCK_HH
#define MCS_CORE_MEMORY_BLOCK_HH


#include "chunk.hh"


namespace mcs
{

  namespace core
  {

    namespace memory
    {


      template<typename Size>
      class block
      {

      public:

        typedef Size size_type;

        typedef chunk<size_type> chunk_type;

        typedef block<size_type> self_type;


      private:

        char* ptr_;


      public:

        block(char* ptr = 0);

        block(char* ptr,
              size_type size);

        char*
        ptr() const;

        size_type
        size() const;

        size_type
        chunk_area() const;

        self_type
        next() const;

        void
        next(self_type next);

        chunk_type
        begin() const;

        chunk_type
        end() const;


      private:

        struct descriptor
        {
          size_type size_;
          char* next_ptr_;

          static size_type& size(self_type block);
          static char*& next_ptr(self_type block);
        };

      };


    }

  }

}


#include "block.cc"
#endif
