
//          Copyright Marc Hofmann 2009 - 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file storage.hh
 *
 * @author Marc Hofmann
 *
 * Work in progress.
 *
 * Derived from Boost Pool.
 */

#ifndef MCS_CORE_MEMORY_STORAGE_HH
#define MCS_CORE_MEMORY_STORAGE_HH


#include "chunk.hh"
#include "block.hh"


namespace mcs
{

  namespace core
  {

    namespace memory
    {


      struct OrderingPolicy
      {

        struct ordered
        {
          static ordered tag;
        };

        struct unordered
        {
          static unordered tag;
        };

      };


      template <typename Policy, typename Size>
      class storage
      {

      public:

        typedef Size size_type;

        typedef chunk<size_type> chunk_type;

        typedef block<size_type> block_type;

        typedef Policy policy_type;


        static policy_type policy_tag;


      private:

        char dummy_chunk_[sizeof(chunk_type)];


      protected:

        size_type chunk_size_;

        chunk_type tail_;


      public:

        storage(size_type chunk_size);


      private:

        storage(const storage&);

        storage&
        operator =(const storage&);


      public:

        bool
        empty() const;

        void
        put_block(block_type block);

        chunk_type
        get_chunk();

        void
        put_chunk(chunk_type chunk);

        chunk_type
        get_chunks(size_type chunk_count);

        void
        put_chunks(chunk_type chunk,
                   size_type chunk_count);

        chunk_type
        begin() const;

        chunk_type
        end() const;


      private:

        void
        put_block_impl(block_type block,
                       OrderingPolicy::unordered);

        void
        put_block_impl(block_type block,
                       OrderingPolicy::ordered);

        void
        put_chunk_impl(chunk_type chunk,
                       OrderingPolicy::unordered);

        void
        put_chunk_impl(chunk_type chunk,
                       OrderingPolicy::ordered);

        chunk_type
        get_chunks_impl(size_type chunk_count,
                        OrderingPolicy::unordered);

        chunk_type
        get_chunks_impl(size_type chunk_count,
                        OrderingPolicy::ordered);

        void
        put_chunks_impl(chunk_type chunk,
                        size_type chunk_count,
                        OrderingPolicy::unordered);

        void
        put_chunks_impl(chunk_type chunk,
                        size_type chunk_count,
                        OrderingPolicy::ordered);

        chunk_type
        make_chunks(chunk_type begin,
                    chunk_type end,
                    chunk_type tail = 0);

        chunk_type
        find_prev(chunk_type chunk);

      };


    }

  }

}


#include "storage.cc"
#endif
