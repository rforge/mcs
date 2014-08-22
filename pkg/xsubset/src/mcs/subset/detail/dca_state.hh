#ifndef MCS_SUBSET_DETAIL_DCA_STATE_HH
#define MCS_SUBSET_DETAIL_DCA_STATE_HH


#include <vector>



namespace mcs    {
namespace subset {
namespace detail {


  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  class DcaTable;


  namespace Preorder {

    template<typename TReal>
    class None;

    template<typename TReal>
    class Full;

    template<typename TReal>
    class Single;

  }



  template<typename TReal>
  class DcaState
  {

    friend Preorder::None<TReal>;
    friend Preorder::Full<TReal>;
    friend Preorder::Single<TReal>;


  private:

    struct Node;


  private:

    std::vector<int>   sWork_;
    std::vector<TReal> rzWork_;
    int                rzLdim_;

    std::vector<Node> nodeStack_;
    Node*             currentNode_;
    Node*             nextNode_;

    TReal currentRss_;

    int nodeCount_;


  private:

    DcaState(int size, int mark, int ldim);


  public:

    DcaState(int size, int mark,
             const int* v, const TReal* rz, int ldrz);

    DcaState(int m, int size, int mark,
             const int* v, const TReal* ay, int lday);

    template<
      template<typename R>
      class TPreorder
      >
    DcaState(int size, int mark,
             const int* v, const TReal* rz, int ldrz,
             const TPreorder<TReal>& p);

    template<
      template<typename R>
      class TPreorder
      >
    DcaState(int m, int size, int mark,
             const int* v, const TReal* ay, int lday,
             const TPreorder<TReal>& p);

    void
    dropColumn(int mark);

    void
    nextNode();

    template<
      template<typename R>
      class TPreorder
      >
    void
    nextNode(const TPreorder<TReal>& p);

    bool
    isDone() const;

    int
    currentSize() const;

    int
    currentMark() const;

    TReal
    minBound() const;

    template<
      template<typename R>
      class TCrit
    >
    TReal
    minBound(int mark, const TCrit<TReal>& c) const;

    template<
      template<typename R>
      class TCrit
    >
    void
    reportSubleading(DcaTable<TReal, TCrit>& table) const;

    int
    nodeCount() const;

  };


}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/dca_state.cc"
#endif
