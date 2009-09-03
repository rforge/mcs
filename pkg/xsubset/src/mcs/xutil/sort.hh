#ifndef _MCS_XUTIL_SORT_HH_
#define _MCS_XUTIL_SORT_HH_


#include <iostream>


#include <functional>
#include <algorithm>


namespace MCS
{

  namespace xutil
  {


    class seq_gen
    {
    private:
      int _cur;

    public:
      seq_gen(int first) : _cur(first)
      {
      }

      int operator ()()
      {
        return _cur++;
      }
    };


    template<typename RandomAccessIter, typename StrictWeakOrdering>
    class index_comp_ftor : public std::binary_function<int, int, bool>
    {
    private:
      RandomAccessIter _first;
      StrictWeakOrdering _comp;

    public:
      index_comp_ftor(RandomAccessIter first, StrictWeakOrdering comp) :
        _first(first), _comp(comp)
      {
      }

      bool operator ()(int i, int j)
      {
        return _comp(_first[i], _first[j]);
      }
    };


    template<typename RandomAccessIter, typename StrictWeakOrdering>
    index_comp_ftor<RandomAccessIter, StrictWeakOrdering>
    index_comp(RandomAccessIter first, StrictWeakOrdering comp)
    {
      return index_comp_ftor<RandomAccessIter, StrictWeakOrdering>(first, comp);
    }


    template<typename RandomAccessIter1, typename RandomAccessIter2>
    void
    permute(RandomAccessIter1 first, int count,
            RandomAccessIter2 iperm)
    {
      for(int i = 0; i < count; ++i)
        {
          int j = iperm[i];
          while (j > i)
            {
              j = iperm[j];
            }
          if (j == i)
            {
              j = iperm[i];
              while (j > i)
                {
                  std::swap(first[i], first[j]);
                  j = iperm[j];
                }
            }
        }
    }


  }

}


#endif
