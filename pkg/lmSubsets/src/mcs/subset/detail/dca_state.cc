#ifndef MCS_SUBSET_DETAIL_DCA_STATE_CC
#define MCS_SUBSET_DETAIL_DCA_STATE_CC


#include "mcs/subset/detail/dca_state.hh"

#include <algorithm>  // std::copy_n
#include <utility>  // std::swap



namespace mcs    {
namespace subset {
namespace detail {


  template<typename TReal>
  class Math;

  template<typename TReal>
  class Lapack;

  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  class DcaTable;



  template<typename TReal>
  struct DcaState<TReal>::Node
  {

    int n;
    int k;

    int*   s;
    TReal* rz;

  };


  template<typename TReal>
  DcaState<TReal>::DcaState(const int size, const int mark, const int ldim) :
    mark_       {mark}                                 ,
    sWork_      ((size - mark + 1) *  size            ),
    rzWork_     ((size - mark + 1) * (size + 1) * ldim),
    rzLdim_     {ldim}                                 ,
    nodeStack_  (size)                                 ,
    currentNode_(nodeStack_.data())                    ,
    nextNode_   (nodeStack_.data())                    ,
    currentRss_ {Math<TReal>::MAX}                     ,
    nodeCount_  {0}
  {
    int*   s  = sWork_.data();
    TReal* rz = rzWork_.data();

    currentNode_->s  = s;
    currentNode_->rz = rz;

    for (int i = 0; i < size; ++i)
      {
        s  += size;
        rz += (size + 1) * rzLdim_;

	nodeStack_[i].s  = s;
        nodeStack_[i].rz = rz;
      }
  }


  template<typename TReal>
  DcaState<TReal>::DcaState(const int size, const int mark, const int* const v,
                            const TReal* const rz, const int ldrz) :
    DcaState(size, mark, size + 1)
  {
    const int n = size;
    const int k = mark;

    ++nextNode_;

    Node* const nxt = nextNode_;

    nxt->n = n;
    nxt->k = k;

    std::copy_n(v, n, nxt->s);
    Lapack<TReal>::lacpy(n + 1, n + 1, rz, ldrz, nxt->rz, rzLdim_);
  }


  template<typename TReal>
  DcaState<TReal>::DcaState(const int m, const int size, const int mark,
                            const int* const v, const TReal* const ay,
                            const int lday) :
    DcaState(size, mark, m)
  {
    const int n = size;
    const int k = mark;

    ++nextNode_;

    Node* const nxt = nextNode_;

    nxt->n = n;
    nxt->k = k;

    std::copy_n(v, n, nxt->s);
    Lapack<TReal>::lacpy(m, n + 1, ay , lday, nxt->rz, rzLdim_);
    Lapack<TReal>::geqrf(m, n + 1, nxt->rz, rzLdim_);
  }


  template<typename TReal>
  template<
    template<typename R>
    class TPreorder
  >
  DcaState<TReal>::DcaState(const int size, const int mark, const int* const v,
                            const TReal* const rz, const int ldrz,
                            const TPreorder<TReal>& preo) :
    DcaState(size, mark, size + 1)
  {
    const int n = size;
    const int k = mark;

    ++nextNode_;

    Node* const tmp = currentNode_;
    Node* const nxt = nextNode_;

    tmp->n = n;
    tmp->k = k;

    std::copy_n(v, n, tmp->s);
    Lapack<TReal>::lacpy(n + 1, n + 1, rz, ldrz, tmp->rz, rzLdim_);

    preo.apply(*tmp, *nxt, rzLdim_);
  }


  template<typename TReal>
  template<
    template<typename R>
    class TPreorder
  >
  DcaState<TReal>::DcaState(const int m, const int size, const int mark,
			    const int* const v, const TReal* const ay,
                            const int lday, const TPreorder<TReal>& preo) :
    DcaState(size, mark, m)
  {
    const int n = size;
    const int k = mark;

    ++nextNode_;

    Node* const tmp = currentNode_;
    Node* const nxt = nextNode_;

    tmp->n = n;
    tmp->k = k;

    std::copy_n(v, n, tmp->s);
    Lapack<TReal>::lacpy(m, n + 1, ay , lday, tmp->rz, rzLdim_);
    Lapack<TReal>::geqrf(m, n + 1, tmp->rz, rzLdim_);

    preo.apply(*tmp, *nxt, rzLdim_);
  }


  template<typename TReal>
  void
  DcaState<TReal>::dropColumn(const int mark)
  {
    ++nextNode_;

    const Node* const cur = currentNode_;
          Node* const nxt = nextNode_;

    const int n = cur->n;
    const int k = mark;

    // dim
    nxt->n = n - 1;
    nxt->k = k;

    // subset
    std::copy_n(cur->s + k + 1, n - k - 1, std::copy_n(cur->s, k, nxt->s));

    // QRD
    Qrz<TReal>::drop(n - k, cur->rz + k * rzLdim_ + k, rzLdim_,
                            nxt->rz + k * rzLdim_ + k, rzLdim_);
  }


  template<typename TReal>
  void
  DcaState<TReal>::nextNode()
  {
    std::swap(*nextNode_, *currentNode_);
    --nextNode_;

    const Node* const cur = currentNode_;
    const int n = cur->n;

    currentRss_ = Math<TReal>::sqr(cur->rz[n * rzLdim_ + n]);

    ++nodeCount_;
  }


  template<typename TReal>
  template<
    template<typename R>
    class TPreorder
    >
  void
  DcaState<TReal>::nextNode(const TPreorder<TReal>& preo)
  {
    preo.apply(*nextNode_, *currentNode_, rzLdim_);
    --nextNode_;

    const Node* const cur = currentNode_;
    const int n = cur->n;

    currentRss_ = Math<TReal>::sqr(cur->rz[n * rzLdim_ + n]);

    ++nodeCount_;
  }


  template<typename TReal>
  bool
  DcaState<TReal>::isFinal() const
  {
    return currentNode_ == nextNode_;
  }


  template<typename TReal>
  int
  DcaState<TReal>::currentSize() const
  {
    return currentNode_->n;
  }


  template<typename TReal>
  int
  DcaState<TReal>::currentMark() const
  {
    return currentNode_->k;
  }


  template<typename TReal>
  TReal
  DcaState<TReal>::minBound() const
  {
    return currentRss_;
  }


  template<typename TReal>
  TReal
  DcaState<TReal>::minBound(const int size, const TReal* const tau) const
  {
    return tau[size - 1] * currentRss_;
  }


  template<typename TReal>
  template<
    template<typename R>
    class TCriterion
  >
  TReal
  DcaState<TReal>::minBound(const int size, const TCriterion<TReal>& crit) const
  {
    return crit.value(size, currentRss_);
  }


  template<typename TReal>
  template<
    template<typename R>
    class TCriterion
  >
  void
  DcaState<TReal>::reportSubleading(DcaTable<TReal,TCriterion>& table) const
  {
    const Node* const cur = currentNode_;

    const int n = cur->n;
    const int k = cur->k;

    table.addSubset(n, currentRss_, cur->s);

    const TReal* z   = cur->rz + n * rzLdim_ + n - 1;
          TReal  rss = currentRss_;
    for (int i = n - 1; i > k; --i, --z)
      {
        rss += Math<TReal>::sqr(*z);

        table.addSubset(i, rss, cur->s);
      }
  }


  template<typename TReal>
  int
  DcaState<TReal>::nodeCount() const
  {
    return nodeCount_;
  }


}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif
