#ifndef MCS_SUBSET_DETAIL_DCA_TABLE_HH
#define MCS_SUBSET_DETAIL_DCA_TABLE_HH


#include "mcs/subset/detail/criteria.hh"



namespace mcs    {
namespace subset {
namespace detail {



  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  class DcaTable
  {

  private:

    int size_;
    int nbest_;

    int*   index_;
    TReal* crit_;
    int*   s_;

    TCrit<TReal> c_;


  public:

    DcaTable(int size, int nbest,
             int* index, TReal* crit, int* subset,
	     const TCrit<TReal>& c);

    TReal
    maxBound() const;

    void
    addSubset(int size, TReal rss, const int* s);

    void
    sortSubsets();

  };



  template<typename TReal>
  class DcaTable<TReal, Criteria::None>
  {

  private:

    int size_;
    int nbest_;

    int* index_;
    TReal* rss_;
    int* s_;


  public:

    DcaTable(int size, int nbest,
             int* index, TReal* rss, int* subset);

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
