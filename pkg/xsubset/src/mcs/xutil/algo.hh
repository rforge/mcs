#ifndef _MCS_XUTIL_ALGO_HH_
#define _MCS_XUTIL_ALGO_HH_


#include <algorithm>
#include <functional>


namespace MCS
{

  namespace xutil
  {


    //--------------------------------------------------------------------------
    // copy_n


    template<typename InputIterator, typename Size,
             typename OutputIterator>
    OutputIterator
    copy_n(InputIterator first, Size count,
           OutputIterator result)
    {
      return std::copy(first, first + count, result);
    };



    //--------------------------------------------------------------------------
    // sort_index


    template<typename RandomAccessIter, typename Size,
             typename StrictWeakOrdering>
    class _index_compare : std::binary_function<Size, Size, bool>
    {
    private:
      RandomAccessIter _first;
      StrictWeakOrdering _comp;

    public:
      _index_compare(RandomAccessIter first,
                     StrictWeakOrdering comp) :
        _first(first), _comp(comp)
      {
      }

      bool operator ()(Size i, Size j)
      {
        return _comp(_first[i], _first[j]);
      }
    };


    template<typename RandomAccessIter1, typename Size,
             typename RandomAccessIter2, typename StrictWeakOrdering>
    void
    sort_index(RandomAccessIter1 first, Size count,
               RandomAccessIter2 index,
               StrictWeakOrdering comp)
    {
      std::sort(index, index + count,
                _index_compare<RandomAccessIter1, Size,
                               StrictWeakOrdering>(first, comp));
    }



    //--------------------------------------------------------------------------
    // permute


    template<typename RandomAccessIter1, typename Size,
             typename RandomAccessIter2>
    void
    permute(RandomAccessIter1 first, Size count,
            RandomAccessIter2 index)
    {
      for (int i = 0; i < count; ++i)
        {
          int j = index[i];
          while (j < i)
            {
              j = index[j];
            }
          std::swap(first[i], first[j]);
        }
    }



    //--------------------------------------------------------------------------
    // rank


    template<typename ForwardIter, typename Size,
             typename RandomAccessIter>
    void
    rank(ForwardIter index, Size count, RandomAccessIter irank)
    {
      for (int i = 0; i < count; ++i)
        {
          irank[*index] = i;
          ++index;
        }
    }


  }

}


#endif
