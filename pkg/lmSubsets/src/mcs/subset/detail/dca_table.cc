#ifndef MCS_SUBSET_DETAIL_DCA_TABLE_CC
#define MCS_SUBSET_DETAIL_DCA_TABLE_CC


#include "mcs/subset/detail/dca_table.hh"

#include <algorithm>  // std::copy_n, std::fill, std::make_heap, std::pop_heap, std::push_heap, std::sort_heap, std::for_each
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
                                       int* const sIndex, TReal* const sRss, TReal* const sCrit, int* const sSize, int* const s,
                                       const TCriterion<TReal>& c) :
    size_  {size} ,
    nbest_ {nbest},
    sIndex_{sIndex},
    sRss_  {sRss}  ,
    sCrit_ {sCrit} ,
    sSize_ {sSize} ,
    s_     {s}    ,
    c_     {c}
  {
    std::iota(sIndex_, sIndex_ + nbest_        , 0               );
    std::fill(sRss_  , sRss_   + nbest_        , Math<TReal>::MAX);
    std::fill(sCrit_ , sCrit_  + nbest_        , Math<TReal>::MAX);
    std::fill(sSize_ , sSize_  + nbest_        , -1              );
    std::fill(s_     , s_      + nbest_ * size_, -1              );
  }


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  TReal
  DcaTable<TReal,TCriterion>::maxBound() const
  {
    const int index  = sIndex_[0];

    return sCrit_[index];
  }


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  void
  DcaTable<TReal,TCriterion>::addSubset(const int size, const TReal rss, const int* s)
  {
    const int index = sIndex_[0];

    const TReal crit = c_.value(size, rss);

    if (sCrit_[index] <= crit) {
      return;
    }

    const auto comp = [=] (const int i, const int j) -> bool {
      return sCrit_[i] < sCrit_[j];
    };

    std::pop_heap(sIndex_, sIndex_ + nbest_, comp);

    const int insert = sIndex_[nbest_ - 1];
    sRss_ [insert] = rss;
    sCrit_[insert] = crit;
    sSize_[insert] = size;
    std::copy_n(s, size, s_ + insert * size_);

    std::push_heap(sIndex_, sIndex_ + nbest_, comp);
  }


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  void
  DcaTable<TReal,TCriterion>::sortSubsets()
  {
    const auto comp = [=] (const int i, const int j) -> bool {
      return sCrit_[i] < sCrit_[j];
    };

    std::sort_heap(sIndex_, sIndex_ + nbest_, comp);

    for (int i = 0; i < nbest_; ++i)
      {
	int* const begin = s_ + i * size_;
	int* const end   = begin + sSize_[i];

	std::sort(begin, end);
      }
  }


  template<typename TReal>
  DcaTable<TReal,Criteria::None>::DcaTable(const int size, const int mark, const int nbest,
                                           int* const sIndex, TReal* const sRss, int* const s) :
    size_  {size}  ,
    mark_  {mark}  ,
    nbest_ {nbest} ,
    sIndex_{sIndex},
    sRss_  {sRss}  ,
    s_     {s}
  {
    std::iota(sIndex_, sIndex_ + (size_ - mark) * nbest_        , 0               );
    std::fill(sRss_  , sRss_   + (size_ - mark) * nbest_        , Math<TReal>::MAX);
    std::fill(s_     , s_      + (size_ - mark) * nbest_ * size_, -1              );
  }


  template<typename TReal>
  TReal
  DcaTable<TReal,Criteria::None>::maxBound(const int size) const
  {
    const int offset = (size - mark_ - 1) * nbest_;
    const int index  = sIndex_[offset];

    return sRss_[index];
  }


  template<typename TReal>
  void
  DcaTable<TReal,Criteria::None>::addSubset(const int size, const TReal rss, const int* s)
  {
    const int offset = (size - mark_ - 1) * nbest_;
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

    const int insert = *(end - 1);
    sRss_[insert] = rss;
    std::copy_n(s, size, s_ + insert * size_);

    std::push_heap(begin, end, comp);
  }


  template<typename TReal>
  void
  DcaTable<TReal,Criteria::None>::sortSubsets()
  {
    const auto comp = [=] (const int i, const int j) -> bool {
      return sRss_[i] < sRss_[j];
    };

    for (int i = 0; i < size_ - mark_; ++i)
      {
        const int offset = i * nbest_;

        int* const begin = sIndex_ + offset;
        int* const end   = begin   + nbest_;

        std::sort_heap(begin, end, comp);

	std::for_each(begin, end, [=] (const int j) {
	    int* const begin = s_ + j * size_;
	    int* const end   = begin + mark_ + i + 1;

	    std::sort(begin, end);
	  });
      }
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif
