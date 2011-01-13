/**
 * @file algo.hh
 */
#ifndef MCS_CORE_UTIL_ALGO_HH
#define MCS_CORE_UTIL_ALGO_HH


namespace mcs
{

  namespace core
  {

    namespace util
    {

      namespace algo
      {


        template<typename InputIterator,
                 typename Size,
                 typename OutputIterator>
        OutputIterator
        copy_n(InputIterator first,
               Size n,
               OutputIterator result);


	template<typename ForwardIterator,
		 typename T>
	void
	iota(ForwardIterator first,
             ForwardIterator last,
             T value);


	template<typename ForwardIterator,
		 typename Size,
		 typename T>
	void
	iota_n(ForwardIterator first,
	       Size n,
	       T value);


	template<typename RandomAccessIterator1,
		 typename RandomAccessIterator2,
		 typename StrictWeakOrdering>
	void
	order(RandomAccessIterator1 first,
              RandomAccessIterator1 last,
              RandomAccessIterator2 value,
              StrictWeakOrdering comp);


	template<typename RandomAccessIterator1,
		 typename Size,
		 typename RandomAccessIterator2,
		 typename StrictWeakOrdering>
	void
	order_n(RandomAccessIterator1 first,
                Size n,
                RandomAccessIterator2 value,
                StrictWeakOrdering comp);


	template<typename RandomAccessIterator1,
		 typename RandomAccessIterator2>
	void
	permute(RandomAccessIterator1 first,
		RandomAccessIterator2 last,
		RandomAccessIterator2 index);


	template<typename RandomAccessIterator1,
		 typename Size,
		 typename RandomAccessIterator2>
	void
	permute_n(RandomAccessIterator1 first,
		  Size n,
		  RandomAccessIterator2 index);


      }

    }

  }

}


#include "algo.cc"
#endif
