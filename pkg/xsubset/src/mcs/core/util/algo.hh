/**
 * @file algo.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_UTIL_ALGO_HH
#define MCS_CORE_UTIL_ALGO_HH


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


	template<typename ForwardIterator,
		 typename Size,
		 typename T>
	void
	iota_n(ForwardIterator first,
	       Size n,
	       T value)
	{
	  for (Size i = 0; i < n; ++i)
	    {
	      *first++ = value++;
	    }
	}


	//--------------------------------------------------------------------------
	// sort_index


	namespace detail
	{


	  template<typename RandomAccessIterator,
		   typename StrictWeakOrdering>
	  class index_compare :
	    std::binary_function<typename std::iterator_traits<RandomAccessIterator>::difference_type,
				 typename std::iterator_traits<RandomAccessIterator>::difference_type,
				 bool>
	  {

	  private:

	    typedef RandomAccessIterator Iter;

	    typedef StrictWeakOrdering Comp;

	    typedef typename std::iterator_traits<RandomAccessIterator>::difference_type Diff;


	    Iter value_;

	    Comp comp_;


	  public:

	    index_compare(const Iter value,
			  const Comp comp) :
	      value_(value),
	      comp_(comp)
	    {
	    }

	    bool operator ()(const Diff i,
			     const Diff j)
	    {
	      return comp_(value_[i], value_[j]);
	    }

	  };


	}


	template<typename RandomAccessIterator1,
		 typename RandomAccessIterator2,
		 typename StrictWeakOrdering>
	void
	index_sort(RandomAccessIterator1 first,
		   RandomAccessIterator1 last,
		   RandomAccessIterator2 value,
		   StrictWeakOrdering comp)
	{
	  typedef typename detail::template index_compare<RandomAccessIterator2, StrictWeakOrdering> Comp;

	  std::sort(first, last, Comp(value, comp));
	}


	template<typename RandomAccessIterator1,
		 typename Size,
		 typename RandomAccessIterator2,
		 typename StrictWeakOrdering>
	void
	index_sort_n(RandomAccessIterator1 first,
		     Size n,
		     RandomAccessIterator2 value,
		     StrictWeakOrdering comp)
	{
	  typedef typename detail::template index_compare<RandomAccessIterator2, StrictWeakOrdering> IndexComp;

	  std::sort(first, first + n, IndexComp(value, comp));
	}



	//--------------------------------------------------------------------------
	// permute


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


	template<typename RandomAccessIterator1,
		 typename RandomAccessIterator2>
	void
	permute(RandomAccessIterator1 first,
		RandomAccessIterator2 last,
		RandomAccessIterator2 index)
	{
	  permute_n(first, last - first, index);
	}


	//--------------------------------------------------------------------------
	// rank


	template<typename InputIterator,
	       typename Size,
	       typename RandomAccessIterator>
	void
	rank_n(InputIterator index,
	       Size n,
	       RandomAccessIterator irank)
	{
	  for (int i = 0; i < n; ++i)
	    {
	      irank[*index] = i;
	      ++index;
	    }
	}


      }

    }

  }

}


#endif
