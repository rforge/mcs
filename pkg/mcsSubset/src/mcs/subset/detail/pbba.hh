#ifndef MCS_SUBSET_DETAIL_PBBA_HH
#define MCS_SUBSET_DETAIL_PBBA_HH



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
    class TCriterion,
    template<typename R>
    class TPreorder
  >
  int
  pbba(DcaState<TReal>& state, DcaTable<TReal,TCriterion>& table,
       const TCriterion<TReal>& c, const TPreorder<TReal>& p);


  template<
    typename TReal,
    template<typename R>
    class TPreorder
   >
  int
  pbba(DcaState<TReal>& state, DcaTable<TReal,Criteria::None>& table,
       const TPreorder<TReal>& p);



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/pbba.cc"
#endif
