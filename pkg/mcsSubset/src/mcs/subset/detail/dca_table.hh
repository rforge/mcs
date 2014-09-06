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
    TReal* sCrit_;
    int*   sSize_;
    int*   s_;

    TCriterion<TReal> c_;


  public:

    DcaTable(int size, int nbest,
             int* sIndex, TReal* sRss, TReal* sCrit, int* sSize, int* s,
	     const TCriterion<TReal>& c);

    TReal
    maxBound() const;

    void
    addSubset(int size, TReal rss, const int* s);

    void
    sortSubsets();

  };



  template<typename TReal>
  class DcaTable<TReal,Criteria::None>
  {

  private:

    int size_;
    int mark_;
    int nbest_;

    int*   sIndex_;
    TReal* sRss_;
    int*   s_;


  public:

    DcaTable(int size, int mark, int nbest,
             int* sIndex, TReal* sRss, int* s);

    TReal
    maxBound(int size) const;

    void
    addSubset(int size, TReal rss, const int* s);

    void
    sortSubsets();

  };



}  // end detail
}  // end subset
}  // end mcs



#include "mcs/subset/detail/dca_table.cc"
#endif
