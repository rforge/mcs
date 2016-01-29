#ifndef MCS_SUBSET_DETAIL_DCA_TABLE_HH
#define MCS_SUBSET_DETAIL_DCA_TABLE_HH


#include "mcs/subset/detail/criteria.hh"



namespace mcs    {
namespace subset {
namespace detail {



  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  class DcaTable
  {

  private:

    int size_;
    int nbest_;

    int*   sIndex_;
    TReal* sRss_;
    TReal* sVal_;
    int*   sWhich_;

    TCriterion<TReal> crit_;


  public:

    DcaTable(int size, int nbest, int* sIndex, TReal* sRss, TReal* sVal,
             int* sWhich, const TCriterion<TReal>& crit);

    TReal
    bestValue() const;

    void
    addSubset(int size, TReal rss, const int* s);

    void
    sortIndex(int* sIndex);

  };



  template<typename TReal>
  class DcaTable<TReal,Criteria::None>
  {

  private:

    int size_;
    int nbest_;

    int*   sIndex_;
    TReal* sRss_;
    int*   sWhich_;


  public:

    DcaTable(int size, int nbest, int* sIndex, TReal* sRss, int* sWhich);

    TReal
    bestValue(int size) const;

    void
    addSubset(int size, TReal rss, const int* s);

    void
    sortIndex(int* sIndex);

  };



}  // end detail
}  // end subset
}  // end mcs



#include "mcs/subset/detail/dca_table.cc"
#endif
