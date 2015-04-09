#ifndef MCS_SUBSET_DETAIL_CRITERIA_HH
#define MCS_SUBSET_DETAIL_CRITERIA_HH


namespace mcs    {
namespace subset {
namespace detail {



  namespace Criteria
  {


    template<typename TReal>
    class None
    {

      TReal value(const int size, const TReal rss) const
      {
	return rss;
      }

    };


  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs


#endif
