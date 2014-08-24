/**
 * @file algo.cc
 */
#ifndef MCS_CORE_UTIL_ALGO_CC
#define MCS_CORE_UTIL_ALGO_CC


#include <algorithm>
#include <iterator>


namespace mcs
{

  namespace core
  {

    namespace util
    {

      namespace algo
      {


	namespace detail
	{


	  template<typename RandomAccessIterator,
		   typename StrictWeakOrdering>
	  class index_compare :
	    std::binary_function<
              typename std::iterator_traits<RandomAccessIterator>::difference_type,
              typename std::iterator_traits<RandomAccessIterator>::difference_type,
              bool
            >
	  {

	  private:

	    RandomAccessIterator value_;

	    StrictWeakOrdering comp_;


	  public:

	    index_compare(const RandomAccessIterator value,
			  const StrictWeakOrdering comp) :
	      value_(value),
	      comp_(comp)
	    {
	    }

	    bool operator ()(const typename std::iterator_traits<
                               RandomAccessIterator>::difference_type i,
			     const typename std::iterator_traits<
                               RandomAccessIterator>::difference_type j)
	    {
	      return comp_(value_[i], value_[j]);
	    }

	  };


	}



        template<typename InputIterator,
                 typename Size,
                 typename OutputIterator>
        OutputIterator
        copy_n(InputIterator first,
               Size n,
               OutputIterator result)
        {
          return std::copy(first, first + n, result);
        }


	template<typename ForwardIterator,
		 typename T>
	void
	iota(ForwardIterator first,
             ForwardIterator last,
             T value)
	{
          while (first != last)
	    {
	      *first++ = value++;
	    }
	}


	template<typename ForwardIterator,
		 typename Size,
		 typename T>
	void
	iota_n(ForwardIterator first,
	       Size n,
	       T value)
	{
	  for (; n > 0; --n)
	    {
	      *first++ = value++;
	    }
	}


	template<typename RandomAccessIterator1,
		 typename RandomAccessIterator2,
		 typename StrictWeakOrdering>
	void
	order(RandomAccessIterator1 first,
              RandomAccessIterator1 last,
              RandomAccessIterator2 value,
              StrictWeakOrdering comp)
	{
	  typedef typename detail::template
            index_compare<RandomAccessIterator2, StrictWeakOrdering> Comp;

	  std::sort(first, last, Comp(value, comp));
	}


	template<typename RandomAccessIterator1,
		 typename Size,
		 typename RandomAccessIterator2,
		 typename StrictWeakOrdering>
	void
	order_n(RandomAccessIterator1 first,
                Size n,
                RandomAccessIterator2 value,
                StrictWeakOrdering comp)
	{
          order(first, first + n, value, comp);
	}


	template<typename RandomAccessIterator1,
		 typename RandomAccessIterator2>
	void
	permute(RandomAccessIterator1 first,
		RandomAccessIterator2 last,
		RandomAccessIterator2 index)
	{
	  permute_n(first, last - first, index);
	}


	template<typename RandomAccessIterator1,
		 typename Size,
		 typename RandomAccessIterator2>
	void
	permute_n(RandomAccessIterator1 first,
		  Size n,
		  RandomAccessIterator2 index)
	{
	  for (Size i = 0; i < n; ++i)
	    {
	      Size j = index[i];
	      while (j < i)
		{
		  j = index[j];
		}
	      std::swap(first[i], first[j]);
	    }
	}


      }

    }

  }

}


#endif
