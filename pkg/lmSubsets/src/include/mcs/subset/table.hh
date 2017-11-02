#ifndef MCS_SUBSET_TABLE_HH
#define MCS_SUBSET_TABLE_HH



#include "mcs/subset/detail/dca_result.hh"



namespace mcs    {
namespace subset {



template<typename Scalar>
using table_best = std::vector<detail::dca_result<Scalar>>;

template<typename Scalar>
using table_all = std::vector<std::vector<detail::dca_result<Scalar>>>;



}  // namespace subset
}  // namespace mcs



#endif
