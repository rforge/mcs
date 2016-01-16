#ifndef MCS_SUBSET_DETAIL_BBA_HH
#define MCS_SUBSET_DETAIL_BBA_HH



namespace mcs    {
namespace subset {
namespace detail {



  template<typename TReal>
  class DcaState;

  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  class DcaTable;



  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  int
  bba(DcaState<TReal>& state, DcaTable<TReal,TCriterion>& table,
      const TCriterion<TReal>& crit);


  template<typename TReal>
  int
  bba(DcaState<TReal>& state, DcaTable<TReal,Criteria::None>& table);



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/bba.cc"
#endif
