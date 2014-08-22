#ifndef MCS_SUBSET_DETAIL_DCA_TABLE_CC
#define MCS_SUBSET_DETAIL_DCA_TABLE_CC


#include "mcs/subset/detail/dca_table.hh"

#include "mcs/subset/detail/criteria.hh"

#include <algorithm>  // std::copy, std::fill, std::make_heap, std::pop_heap, std::push_heap, std::sort_heap
#include <numeric>  // std::iota



namespace mcs    {
namespace subset {
namespace detail {


  template<typename TReal>
  class Math;



  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  DcaTable<TReal, TCrit>::DcaTable(const int size, const int nbest,
				   int* const index, TReal* const crit, int* const s,
                                   const TCrit<TReal>& c) :
    size_ {size} ,
    nbest_{nbest},
    index_{index},
    crit_ {crit} ,
    s_    {s}    ,
    c_    {c}
  {
    std::iota(index_, index_ + nbest_        , 0               );
    std::fill(crit_ , crit_  + nbest_        , Math<TReal>::MAX);
    std::fill(s_    , s_     + nbest_ * size_, -1              );
  }


  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  TReal
  DcaTable<TReal, TCrit>::maxBound() const
  {
    return crit_[0];
  }


  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  void
  DcaTable<TReal, TCrit>::addSubset(const int size, const TReal rss, const int* s)
  {
    const TReal crit = c_.value(size, rss);

    if (crit_[0] <= crit) {
      return;
    }

    const auto comp = [&] (const int i, const int j) -> bool {
      return crit_[i] < crit_[j];
    };

    std::pop_heap(index_, index_ + nbest_, comp);

    const int insert = index_[nbest_ - 1];
    crit_[insert] = crit;
    std::copy(s, s + size, s_ + insert * size_);

    std::push_heap(index_, index_ + nbest_, comp);
  }


  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  void
  DcaTable<TReal, TCrit>::sortSubsets()
  {
    const auto comp = [&] (const int i, const int j) -> bool {
      return crit_[i] < crit_[j];
    };

    std::sort(index_, index_ + nbest_, comp);
  }


  template<typename TReal>
  DcaTable<TReal, Criteria::None>::DcaTable(const int size, const int nbest,
					    int* const index, TReal* const rss, int* const s) :
    size_ {size} ,
    nbest_{nbest},
    index_{index},
    rss_  {rss}  ,
    s_    {s}
  {
    std::iota(index_, index_ + size_ * nbest_        , 0               );
    std::fill(rss_  , rss_   + size_ * nbest_        , Math<TReal>::MAX);
    std::fill(s_    , s_     + size_ * nbest_ * size_, -1              );
  }


  template<typename TReal>
  TReal
  DcaTable<TReal, Criteria::None>::maxBound(const int size) const
  {
    const int offset = (size - 1) * nbest_;
    const int index  = index_[offset];

    return rss_[index];
  }


  template<typename TReal>
  void
  DcaTable<TReal, Criteria::None>::addSubset(const int size, const TReal rss, const int* s)
  {
    const int offset = (size - 1) * nbest_;
    const int index  = index_[offset];

    if (rss_[index] <= rss) {
      return;
    }

    const auto comp = [&] (const int i, const int j) -> bool {
      return rss_[i] < rss_[j];
    };

    int* const begin = index_ + offset;
    int* const end   = begin  + nbest_;

    std::pop_heap(begin, end, comp);

    const int insert = *(end - 1);
    rss_[insert] = rss;
    std::copy(s, s + size, s_ + insert * size_);

    std::push_heap(begin, end, comp);
  }


  template<typename TReal>
  void
  DcaTable<TReal, Criteria::None>::sortSubsets()
  {
    const auto comp = [&] (const int i, const int j) -> bool {
      return rss_[i] < rss_[j];
    };

    for (int i = 0; i < size_; ++i)
      {
        const int offset = i * nbest_;

        int* const begin = index_ + offset;
        int* const end   = begin  + nbest_;

        std::sort_heap(begin, end, comp);
      }
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif
