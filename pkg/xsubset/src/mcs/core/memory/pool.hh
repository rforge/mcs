
//          Copyright Marc Hofmann 2009 - 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file pool.hh
 *
 * @author Marc Hofmann
 *
 * Work in progress.
 *
 * Derived from Boost Pool.
 */

#ifndef MCS_CORE_MEMORY_MEMORY_POOL_HH
#define MCS_CORE_MEMORY_MEMORY_POOL_HH


#include "block.hh"
#include "storage.hh"


namespace mcs
{

  namespace core
  {

    namespace memory
    {


      template <typename Allocator>
      class memory_pool:
      {

      public:

        typedef Allocator allocator_type;

        typedef typename allocator_type::size_type size_type;

        typedef typename allocator_type::difference_type difference_type;

        typedef block<size_type> block_type;

        typedef chunk<size_type> chunk_type;

        typedef storage<size_type> storage_type;


      private:

        static const min_block_size;


      protected:

        storage_type store_;

        block_type tail_;


      public:

        explicit
        pool(size_type element_size);

        ~pool();

        bool
        release_orderered_memory();

        bool
        purge_memory();

        char*
        alloc_chunk();

        char*
        alloc_ordered_chunk();

        char*
        alloc_ordered_chunks(size_type chunk_count);

        void
        free_chunk(char* ptr);

        void
        free_ordered_chunk(void* ptr);

        void
        free_chunks(char* ptr,
                    size_type chunk_count);

        void
        free_ordered_chunks(char* ptr,
                            size_type chunk_count);

        chunk_type
        chunks_begin() const;

        chunk_type
        chunks_end() const;

      };



    }

  }

}


#endif
