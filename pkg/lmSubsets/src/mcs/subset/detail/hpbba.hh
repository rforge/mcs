#ifndef MCS_SUBSET_DETAIL_HPBBA_HH
#define MCS_SUBSET_DETAIL_HPBBA_HH



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
  hpbba(DcaState<TReal>& state, DcaTable<TReal,TCriterion>& table,
	const TCriterion<TReal>& crit, const TPreorder<TReal>& preo,
	const TReal tau);


  template<
    typename TReal,
    template<typename R>
    class TPreorder
   >
  int
  hpbba(DcaState<TReal>& state, DcaTable<TReal,Criteria::None>& table,
	const TPreorder<TReal>& preo, const TReal* const tau);




}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/hpbba.cc"
#endif
