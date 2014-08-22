#ifndef MCS_SUBSET_DETAIL_PREORDER_CC
#define MCS_SUBSET_DETAIL_PREORDER_CC


#include "mcs/subset/detail/preorder.hh"

#include <numeric>  // std::iota
#include <algorithm>  // std::sort, std::copy_n, std::max_element
#include <utility>  // std::swap



namespace mcs    {
namespace subset {
namespace detail {



  template<typename TReal>
  int
  Preorder::None<TReal>::pmin() const
  {
    return Math<int>::MAX;
  }


  template<typename TReal>
  void
  Preorder::None<TReal>::apply(typename DcaState<TReal>::Node& in,
                               typename DcaState<TReal>::Node& out,
                               const int rzLdim) const
  {
    std::swap(in, out);
  }



  template<typename TReal>
  Preorder::Full<TReal>::Full(int size, int pmin) :
    pmin_{pmin},
    givens_(size + 1),
    bounds_(size + 1),
    permut_ (size + 1)
  {
  }


  template<typename TReal>
  int
  Preorder::Full<TReal>::pmin() const
  {
    return pmin_;
  }


  template<typename TReal>
  void
  Preorder::Full<TReal>::apply(typename DcaState<TReal>::Node& in,
                               typename DcaState<TReal>::Node& out,
                               const int rzLdim) const
  {
    const int n = in.n;
    const int k = in.k;

    // preorder?
    if (n - k <= pmin_)
      {
        std::swap(in, out);

        return;
      }

    // compute bounds
    TReal* const b = bounds_.data();

    Qrz<TReal>::bounds(n - k, in.rz + k * rzLdim + k, rzLdim, b, givens_);

    // compute permutation
    int* const pi = permut_.data();

    std::iota(pi, pi + n - k, 0);
    std::sort(pi, pi + n - k, [&] (const int i, const int j) {
        return b[i] > b[j];
      });

    // dim
    out.n = n;
    out.k = k;

    // permute subset
    std::copy_n(in.s, k, out.s);
    Algo::permute_n(pi, n - k, in.s + k, out.s + k);

    // permute QRD
    Qrz<TReal>::permute(n - k,  in.rz + k * rzLdim + k, rzLdim , pi,
                               out.rz + k * rzLdim + k, rzLdim);
  }



  template<typename TReal>
  Preorder::Single<TReal>::Single(int size, int pmin) :
    pmin_{pmin},
    givens_(size + 1),
    bounds_(size + 1)
  {
  }


  template<typename TReal>
  int
  Preorder::Single<TReal>::pmin() const
  {
    return pmin_;
  }


  template<typename TReal>
  void
  Preorder::Single<TReal>::apply(typename DcaState<TReal>::Node& in,
                                 typename DcaState<TReal>::Node& out,
                                 const int rzLdim) const
  {
    const int n = in.n;
    const int k = in.k;

    std::swap(in, out);

    // preorder?
    if (n - k <= pmin_)
      {
        return;
      }

    // bounds
    TReal* const b = bounds_.data();

    Qrz<TReal>::bounds(n - k, out.rz + k * rzLdim + k, rzLdim, b, givens_);

    // find max
    const auto mx = std::max_element(b, b + n - k);
    const int j = mx - b;

    // preorder subset
    std::swap(out.s[k], out.s[k + j]);

    // preorder QRD
    Qrz<TReal>::permute1(n - k, j, out.rz + k * rzLdim + k, rzLdim);
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif
