#ifndef MCS_SUBSET_DETAIL_DCA_TABLE_CC
#define MCS_SUBSET_DETAIL_DCA_TABLE_CC


#include "mcs/subset/detail/dca_table.hh"

#include <algorithm>  // std::copy_n, std::fill, std::make_heap, std::pop_heap, std::push_heap, std::sort, std::for_each
#include <numeric>  // std::iota


namespace mcs    {
namespace subset {
namespace detail {


  template<typename TReal>
  class Math;



  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  DcaTable<TReal,TCriterion>::DcaTable(const int size, const int nbest,
                                       int* const sIndex, TReal* const sRss,
                                       TReal* const sVal, int* const sWhich,
                                       const TCriterion<TReal>& crit) :
    size_  {size}  ,
    nbest_ {nbest} ,
    sIndex_{sIndex},
    sRss_  {sRss}  ,
    sVal_  {sVal}  ,
    sWhich_{sWhich},
    crit_  {crit}
  {
    std::iota(sIndex_, sIndex_ + nbest_              , 0               );
    std::fill(sRss_  , sRss_   + nbest_              , Math<TReal>::MAX);
    std::fill(sVal_  , sVal_   + nbest_              , Math<TReal>::MAX);
    std::fill(sWhich_, sWhich_ + nbest_ * (size_ + 1), -1              );
  }


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  TReal
  DcaTable<TReal,TCriterion>::bestValue() const
  {
    const int index  = sIndex_[0];

    return sVal_[index];
  }


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  void
  DcaTable<TReal,TCriterion>::addSubset(const int size, const TReal rss,
                                        const int* s)
  {
    const int index = sIndex_[0];

    const TReal val = crit_.value(size, rss);

    if (sVal_[index] <= val) {
      return;
    }

    const auto comp = [=] (const int i, const int j) -> bool {
      return sVal_[i] < sVal_[j];
    };

    std::pop_heap(sIndex_, sIndex_ + nbest_, comp);

    int insert = sIndex_[nbest_ - 1];
    sRss_ [insert] = rss;
    sVal_ [insert] = val;

    insert *= (size_ + 1);
    sWhich_[insert] = size;
    std::copy_n(s, size, sWhich_ + insert + 1);

    std::push_heap(sIndex_, sIndex_ + nbest_, comp);
  }


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  void
  DcaTable<TReal,TCriterion>::sortIndex(int* const sIndex)
  {
    const auto comp = [=] (const int i, const int j) -> bool {
      return sVal_[i] < sVal_[j];
    };

    std::sort(sIndex, sIndex + nbest_, comp);
  }


  template<typename TReal>
  DcaTable<TReal,Criteria::None>::DcaTable(const int size, const int nbest,
                                           int* const sIndex, TReal* const sRss,
                                           int* const sWhich) :
    size_  {size}  ,
    nbest_ {nbest} ,
    sIndex_{sIndex},
    sRss_  {sRss}  ,
    sWhich_{sWhich}
  {
    std::iota(sIndex_, sIndex_ + size_ * nbest_              , 0               );
    std::fill(sRss_  , sRss_   + size_ * nbest_              , Math<TReal>::MAX);
    std::fill(sWhich_, sWhich_ + size_ * nbest_ * (size_ + 1), -1              );
  }


  template<typename TReal>
  TReal
  DcaTable<TReal,Criteria::None>::bestValue(const int size) const
  {
    const int offset = (size - 1) * nbest_;
    const int index  = sIndex_[offset];

    return sRss_[index];
  }


  template<typename TReal>
  void
  DcaTable<TReal,Criteria::None>::addSubset(const int size, const TReal rss,
                                            const int* s)
  {
    const int offset = (size - 1) * nbest_;
    const int index  = sIndex_[offset];

    if (sRss_[index] <= rss) {
      return;
    }

    const auto comp = [=] (const int i, const int j) -> bool {
      return sRss_[i] < sRss_[j];
    };

    int* const begin = sIndex_ + offset;
    int* const end   = begin   + nbest_;

    std::pop_heap(begin, end, comp);

    int insert = *(end - 1);
    sRss_[insert] = rss;

    insert *= (size_ + 1);
    sWhich_[insert] = size;
    std::copy_n(s, size, sWhich_ + insert + 1);

    std::push_heap(begin, end, comp);
  }


  template<typename TReal>
  void
  DcaTable<TReal,Criteria::None>::sortIndex(int* const sIndex)
  {
    const auto comp = [=] (const int i, const int j) -> bool {
      return sRss_[i] < sRss_[j];
    };

    for (int i = 0; i < size_; ++i)
      {
        const int offset = i * nbest_;

        int* const begin = sIndex_ + offset;
        int* const end   = begin   + nbest_;

        std::sort(begin, end, comp);
      }
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif
