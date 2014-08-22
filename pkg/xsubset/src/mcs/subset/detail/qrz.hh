#ifndef MCS_SUBSET_DETAIL_QRZ_HH
#define MCS_SUBSET_DETAIL_QRZ_HH


#include <vector>



namespace mcs    {
namespace subset {
namespace detail {


  template<typename TReal>
  class Givens;



  template<typename TReal>
  struct Qrz
  {


    static
    void
    drop(int n, const TReal* r , int ldr ,
                      TReal* rr, int ldrr);


    template<typename OutputIterator>
    static
    void
    bounds(int n, const TReal* rz, int ldrz,
           OutputIterator out, std::vector<Givens<TReal>>& work);


    template<typename InputIterator>
    static
    void
    permute(int n, const TReal* rz, int ldrz , InputIterator pi,
                         TReal* sz, int ldsz);


    static
    void
    permute2(int n, const TReal* rz, int ldrz , int j,
                          TReal* sz, int ldsz);


    static
    void
    permute1(int n, int j, TReal* rz, int ldrz);


    static
    void
    permute2(int n, int j, TReal* rz, int ldrz);


  };



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/qrz.cc"
#endif
