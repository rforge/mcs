#ifndef MCS_SUBSET_DETAIL_DCA_HH
#define MCS_SUBSET_DETAIL_DCA_HH



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
  dca(DcaState<TReal>& state, DcaTable<TReal,TCriterion>& table);




}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/dca.cc"
#endif
