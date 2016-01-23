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
<<<<<<< .mine
       const TCriterion<TReal>& crit, const TReal tau);
=======
       const TCriterion<TReal>& crit, const TPreorder<TReal>& preo,
       const TReal tau);
>>>>>>> .r101


  template<
    typename TReal
   >
  int
  hbba(DcaState<TReal>& state, DcaTable<TReal,Criteria::None>& table,
<<<<<<< .mine
       const TReal* const tau);
=======
       const TPreorder<TReal>& preo, const TReal* const tau);
>>>>>>> .r101




}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/hbba.cc"
#endif
