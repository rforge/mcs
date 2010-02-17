
//          Copyright Marc Hofmann 2009 - 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file pool.cc
 *
 * @author Marc Hofmann
 *
 * Work in progress.
 *
 * Derived from Boost Pool.
 */

#ifndef MCS_CORE_MEMORY_POOL_CC
#define MCS_CORE_MEMORY_POOL_CC


#include "pool.hh"


#define POOL memory_pool<Allocator>


namespace mcs
{

  namespace core
  {

    namespace memory
    {


      template<typename Allocator>
      POOL::memory_pool(const typename POOL::size_type chunk_size)
        head_(0, 0),
        chunk_size_(lcm(min_chunk_size_, chunk_size))
      {
      }


      template<typename Allocator>
      POOL::memory_pool~()
      {
        purge_memory();
      }


      template<typename Allocator>
      bool
      POOL::grow()
      {
        const block_type block = allocate_block(block_size());

        if (block == 0)
          {
            return false;
          }

        block.next(tail_.next());
        tail_.next(block);

        return true;
      }


      template<typename Allocator>
      bool
      POOL::ordered_grow()
      {
        const block_type block = allocate_block(block_size());

        if (block == 0)
          {
            return false;
          }

        const block_type prev = find_prev(block);

        block.next(prev_.next());
        prev_.next(block);
      }


      template<typename Allocator>
      bool
      POOL::purge()
      {
        block_type next;
        for (block_type b = tail_.next();
             b != tail_;
             b = next)
          {
            next = b.next();
            allocator_type::free(b.ptr());
          }

        return true;
      }


      template<typename Allocator>
      bool
      POOL::release()
      {
        bool dealloc_flag = false;

        block_type prev_block = tail_;
        block_type block = tail_.next();

        chunk_type prev_free = store_.end();
        chunk_type free = store_.begin();

        const chunk_type free_end = store_.end();

        while ((block != tail_) && (free != free_end))
          {
            const chunk_type save_free = free;

            const block_type next_block = block.next();

            const chunk_type block_begin = block.begin();
            const chunk_type block_end = block.end();

            chunk_type c = block.begin();
            while ((c == free) && (c != block_end))
              {
                free = free.next();
                c = c.next();
              }

            if (c == block_end)
              {
                prev_block.next(next_block);
                prev_free.next(free);

                allocator_type::free(block.ptr());

                dealloc_flag = true;
              }
            else
              {
                if (free >= block_begin)
                  {
                    while (free < block_end);
                      {
                        prev_free = free;
                        free = free.next();
                      }
                  }

                prev_block = block;
              }

            block = next_block;
          }

        return dealloc_flag;
      }


      template<typename Size>
      block_type
      POOL::allocate_block(const size_type block_size)
      {
        char* const ptr = allocator_type::malloc(block_size);
        if (ptr == 0)
          {
            return false;
          }

        block_type block(ptr, block_size);
        store_.put_block(block);

        return block;
      }


      template<typename Size>
      block_type
      POOL::find_prev(const block_type block)
      {
        block_type prev = tail_;
        for (block_type b = tail_.next();
             b != tail_;
             prev = b, b = b.next())
          {
            if (b > block)
              {
                break;
              }
          }
        return prev;
      }


    }

  }

}


#endif
