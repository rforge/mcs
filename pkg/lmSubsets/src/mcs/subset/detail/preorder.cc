#ifndef MCS_SUBSET_DETAIL_PREORDER_CC
#define MCS_SUBSET_DETAIL_PREORDER_CC


#include "mcs/subset/detail/preorder.hh"

#include <numeric>  // std::iota
#include <algorithm>  // std::sort, std::copy_n, std::max_element, std::rotate
#include <utility>  // std::swap
#include <iterator>  // std::reverse_iterator



namespace mcs    {
namespace subset {
namespace detail {



  template<typename TReal>
  void
  Preorder::None<TReal>::apply(typename DcaState<TReal>::Node& in,
                               typename DcaState<TReal>::Node& out,
                               const int rzLdim) const
  {
    std::swap(in, out);
  }



  template<typename TReal>
  Preorder::Full<TReal>::Full(int size) :
    Full(size, -1)
  {
  }


  template<typename TReal>
  Preorder::Full<TReal>::Full(int size, int pmin) :
    Full(size, pmin, -1)
  {
  }


  template<typename TReal>
  Preorder::Full<TReal>::Full(int size, int pmin, int pmax) :
    pmin_{pmin},
    pmax_{pmax},
    givens_(size + 1),
    bounds_(size + 1),
    permut_ (size + 1)
  {
  }


  template<typename TReal>
  void
  Preorder::Full<TReal>::apply(typename DcaState<TReal>::Node& in,
                               typename DcaState<TReal>::Node& out,
                               const int rzLdim) const
  {
    const int n = in.n;
    const int k = in.k;
    const int p = n - k;

    // preorder?
    if (((pmin_ >= 0) && (p < pmin_)) ||
        ((pmax_ >= 0) && (p > pmax_)))
      {
        std::swap(in, out);

        return;
      }

    // compute bounds
    TReal* const b = bounds_.data();

    Qrz<TReal>::bounds(p, in.rz + k * rzLdim + k, rzLdim, b, givens_);

    // compute permutation
    int* const pi = permut_.data();

    std::iota(pi, pi + p, 0);
    std::sort(pi, pi + p, [&] (const int i, const int j) {
        return b[i] > b[j];
      });

    // dim
    out.n = n;
    out.k = k;

    // permute subset
    std::copy_n(in.s, k, out.s);
    Algo::permute_n(pi, p, in.s + k, out.s + k);

    // permute QRD
    Qrz<TReal>::permute(p,  in.rz + k * rzLdim + k, rzLdim , pi,
                           out.rz + k * rzLdim + k, rzLdim);
  }



  template<typename TReal>
  Preorder::Single1<TReal>::Single1(int size) :
    Single1(size, -1)
  {
  }


  template<typename TReal>
  Preorder::Single1<TReal>::Single1(int size, int pmin) :
    Single1(size, pmin, -1)
  {
  }


  template<typename TReal>
  Preorder::Single1<TReal>::Single1(int size, int pmin, int pmax) :
    pmin_{pmin},
    pmax_{pmax},
    givens_(size + 1),
    bounds_(size + 1)
  {
  }


  template<typename TReal>
  void
  Preorder::Single1<TReal>::apply(typename DcaState<TReal>::Node& in,
                                  typename DcaState<TReal>::Node& out,
                                  const int rzLdim) const
  {
    const int n = in.n;
    const int k = in.k;
    const int p = n - k;

    std::swap(in, out);

    // preorder?
    if (((pmin_ >= 0) && (p < pmin_)) ||
        ((pmax_ >= 0) && (p > pmax_)))
      {
        return;
      }

    // bounds
    TReal* const b = bounds_.data();

    Qrz<TReal>::bounds(p, out.rz + k * rzLdim + k, rzLdim, b, givens_);

    // find max
    const auto mx = std::max_element(b, b + p);
    const int j = mx - b;

    // preorder v1
    std::swap(out.s[k], out.s[k + j]);
    Qrz<TReal>::permute1(p, j, out.rz + k * rzLdim + k, rzLdim);
  }



  template<typename TReal>
  Preorder::Single2<TReal>::Single2(int size) :
    Single2(size, -1)
  {
  }


  template<typename TReal>
  Preorder::Single2<TReal>::Single2(int size, int pmin) :
    Single2(size, pmin, -1)
  {
  }


  template<typename TReal>
  Preorder::Single2<TReal>::Single2(int size, int pmin, int pmax) :
    pmin_{pmin},
    pmax_{pmax},
    givens_(size + 1),
    bounds_(size + 1)
  {
  }


  template<typename TReal>
  void
  Preorder::Single2<TReal>::apply(typename DcaState<TReal>::Node& in,
                                  typename DcaState<TReal>::Node& out,
                                  const int rzLdim) const
  {
    const int n = in.n;
    const int k = in.k;
    const int p = n - k;

    std::swap(in, out);

    // preorder?
    if (((pmin_ >= 0) && (p < pmin_)) ||
        ((pmax_ >= 0) && (p > pmax_)))
      {
        return;
      }

    // bounds
    TReal* const b = bounds_.data();

    Qrz<TReal>::bounds(p, out.rz + k * rzLdim + k, rzLdim, b, givens_);

    // find max
    const auto mx = std::max_element(b, b + p);
    const int j = mx - b;

    // preorder v2
    std::rotate(std::reverse_iterator<int*>(out.s + k + j + 1),
                std::reverse_iterator<int*>(out.s + k + j),
                std::reverse_iterator<int*>(out.s + k));
    Qrz<TReal>::permute2(p, j, out.rz + k * rzLdim + k, rzLdim);
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs


#endif
