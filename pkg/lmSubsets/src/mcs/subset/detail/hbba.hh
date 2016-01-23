#ifndef MCS_SUBSET_DETAIL_HBBA_HH
#define MCS_SUBSET_DETAIL_HBBA_HH



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
  hbba(DcaState<TReal>& state, DcaTable<TReal,TCriterion>& table,
       const TCriterion<TReal>& crit, const TReal tau);


  template<
    typename TReal
   >
  int
  hbba(DcaState<TReal>& state, DcaTable<TReal,Criteria::None>& table,
       const TReal* const tau);




}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/hbba.cc"
#endif
