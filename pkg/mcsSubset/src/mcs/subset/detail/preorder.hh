#ifndef MCS_SUBSET_DETAIL_PREORDER_HH
#define MCS_SUBSET_DETAIL_PREORDER_HH


#include <vector>



namespace mcs    {
namespace subset {
namespace detail {


  template<typename TReal>
  class DcaState;



  namespace Preorder
  {


    template<typename TReal>
    class None
    {

    public:

      void
      apply(typename DcaState<TReal>::Node& in,
            typename DcaState<TReal>::Node& out,
            int rzLdim) const;

    };



    template<typename TReal>
    class Full
    {

    private:

      int pmin_;

      mutable std::vector<Givens<TReal>> givens_;
      mutable std::vector<TReal>         bounds_;
      mutable std::vector<int>           permut_;


    public:

      Full(int size, int pmin);

      void
      apply(typename DcaState<TReal>::Node& in,
            typename DcaState<TReal>::Node& out,
            int rzLdim) const;

    };



    template<typename TReal>
    class Single
    {

    private:

      int pmin_;

      mutable std::vector<Givens<TReal>> givens_;
      mutable std::vector<TReal>         bounds_;


    public:

      Single(int size, int pmin);

      void
      apply(typename DcaState<TReal>::Node& in,
            typename DcaState<TReal>::Node& out,
            int rzLdim) const;

    };


  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs


#include "mcs/subset/detail/preorder.cc"
#endif
