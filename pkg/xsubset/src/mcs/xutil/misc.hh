#ifndef _MCS_XUTIL_MISC_HH_
#define _MCS_XUTIL_MISC_HH_


#include <algorithm>


namespace MCS
{

  namespace xutil
  {


    template<typename InputIterator, typename Size,
             typename OutputIterator>
    OutputIterator
    copy_n(InputIterator first, Size count,
           OutputIterator result)
    {
      return std::copy(first, first + count, result);
    };


    template<typename RandomAccessIter, typename StrictWeakOrdering>
    class __index_compare : std::binary_function<int, int, bool>
    {
    private:
      RandomAccessIter _first;
      StrictWeakOrdering _comp;

    public:
      index_compare(RandomAccessIter first,
                    StrictWeakOrdering comp) :
        _first(first), _comp(comp)
      {
      }

      bool operator ()(int i, int j)
      {
        return _comp(_first[i], _first[j]);
      }
    };



    template<typename RandomAccessIter1, typename RandomAccessIter2,
             typename StrictWeakOrdering>
    void
    sort_index(RandomAccessIter1 first, int count,
               RandomAccessIter2 iperm,
               StrictWeakOrdering comp)
    {
      std::sort(iperm, iperm + count,
                __index_compare<RandomAccessIter1,StrictWeakOrdering>(first, comp));
    }


  }

}


#endif
