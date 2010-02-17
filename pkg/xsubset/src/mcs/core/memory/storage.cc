
//          Copyright Marc Hofmann 2009 - 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file storage.cc
 *
 * @author Marc Hofmann
 *
 * Work in progress.
 *
 * Derived from Boost Pool.
 */

#ifndef MCS_CORE_MEMORY_STORAGE_CC
#define MCS_CORE_MEMORY_STORAGE_CC


#include "storage.hh"


#define STORAGE storage<Policy, Size>


namespace mcs
{

  namespace core
  {

    namespace memory
    {


      template<typename Policy, typename Size>
      STORAGE::storage(typename STORAGE::size_type chunk_size) :
        tail_(0),
        chunk_size_(chunk_size)
      {
        tail_ = dummy_chunk_;
        tail_.next(tail_);
      }


      template<typename Policy, typename Size>
      bool
      STORAGE::empty() const
      {
        return tail_.next() == tail_;
      }


      template<typename Policy, typename Size>
      void
      STORAGE::put_block(const typename STORAGE::block_type block)
      {
        put_block_impl(block, policy_tag);
      }



      template<typename Policy, typename Size>
      typename STORAGE::chunk_type
      STORAGE::get_chunk()
      {
        const chunk_type chunk = tail_.next();
        tail_.next(chunk.next());
        return chunk;
      }


      template<typename Policy, typename Size>
      void
      STORAGE::put_chunk(const typename STORAGE::chunk_type chunk)
      {
        put_chunk_impl(chunk, policy_tag);
      }


      template<typename Policy, typename Size>
      typename STORAGE::chunk_type
      STORAGE::get_chunks(const typename STORAGE::size_type chunk_count)
      {
        get_chunks_impl(chunk_count, policy_tag);
      }


      template<typename Policy, typename Size>
      void
      STORAGE::put_chunks(const typename STORAGE::chunk_type chunk,
                          const typename STORAGE::size_type chunk_count)
      {
        put_chunks_impl(chunk, chunk_count, policy_tag);
      }


      template<typename Policy, typename Size>
      void
      STORAGE::put_block_impl(const typename STORAGE::block_type block,
                              const OrderingPolicy::unordered)
      {
        tail_.next(make_chunks(block.begin(), block.end(), tail_.next()));
      }


      template<typename Policy, typename Size>
      void
      STORAGE::put_block_impl(const typename STORAGE::block_type block,
                              const OrderingPolicy::ordered)
      {
        chunk_type prev = find_prev(block.begin());
        prev.next(make_chunks(block.begin(), block.end(), prev.next()));
      }


      template<typename Policy, typename Size>
      void
      STORAGE::put_chunk_impl(const typename STORAGE::chunk_type chunk,
                              const OrderingPolicy::unordered)
      {
        chunk.next(tail_.next());
        tail_.next(chunk);
      }


      template<typename Policy, typename Size>
      void
      STORAGE::put_chunk_impl(const typename STORAGE::chunk_type chunk,
                              const OrderingPolicy::ordered)
      {
        chunk_type prev = find_prev(chunk);
        chunk.next(prev.next());
        prev.next(chunk);
      }


      template<typename Policy, typename Size>
      typename STORAGE::chunk_type
      STORAGE::get_chunks_impl(const typename STORAGE::size_type chunk_count,
                               const OrderingPolicy::unordered)
      {
        return 0;
      }


      template<typename Policy, typename Size>
      typename STORAGE::chunk_type
      STORAGE::get_chunks_impl(const typename STORAGE::size_type chunk_count,
                               const OrderingPolicy::ordered)
      {
        if(chunk_count == 0)
          {
            return 0;
          }

        chunk_type prev = tail_;
        chunk_type first = tail_.next();
        chunk_type last;

        do
          {
            if (first == tail_)
              {
                return 0;
              }

            last = first;

            for (int n = chunk_count; n != 1; --n)
              {
                const chunk_type next = last.next();

                if ((next - last) != chunk_size_)
                  {
                    prev = last;
                    first = next;
                    last = 0;
                    break;
                  }

                last = next;
              }
          }
        while (last == 0);

        prev.next(last.next());

        return first;
      }


      template<typename Policy, typename Size>
      void
      STORAGE::put_chunks_impl(const typename STORAGE::chunk_type chunk,
                               const typename STORAGE::size_type chunk_count,
                               const OrderingPolicy::unordered policy_tag)
      {
        if(chunk_count != 0)
          {
            put_block_impl(block_type(chunk.ptr(),
                                      chunk_count * chunk_size_),
                           policy_tag);
          }
      }


      template<typename Policy, typename Size>
      void
      STORAGE::put_chunks_impl(const typename STORAGE::chunk_type chunk,
                               const typename STORAGE::size_type chunk_count,
                               const OrderingPolicy::ordered policy_tag)
      {
        if(chunk_count != 0)
          {
            put_block_impl(block_type(chunk.ptr(),
                                      chunk_count * chunk_size_),
                           policy_tag);
          }
      }


      template<typename Policy, typename Size>
      typename STORAGE::chunk_type
      STORAGE::make_chunks(const typename STORAGE::chunk_type begin,
                           const typename STORAGE::chunk_type end,
                           const typename STORAGE::chunk_type tail)
      {
        chunk_type prev = begin;
        for (chunk_type c = begin + chunk_size_;
             c < end;
             c += chunk_size_)
          {
            prev.next(c);
            prev = c;
          }
        prev.next(tail);

        return begin;
      }


      template<typename Policy, typename Size>
      typename STORAGE::chunk_type
      STORAGE::find_prev(const typename STORAGE::chunk_type chunk)
      {
        chunk_type prev = tail_;
        for (chunk_type c = tail_.next();
             c != tail_;
             prev = c, c = c.next())
          {
            if (c > chunk)
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
