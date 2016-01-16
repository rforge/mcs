#ifndef MCS_SUBSET_DETAIL_DCA_STATE_HH
#define MCS_SUBSET_DETAIL_DCA_STATE_HH


#include <vector>



namespace mcs    {
namespace subset {
namespace detail {


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  class DcaTable;



  template<typename TReal>
  class DcaState
  {

  public:

    struct Node;


  private:

    int mark_;

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

    DcaState(int size, int mark, const int* v,
             const TReal* rz, int ldrz);

    DcaState(int m, int size, int mark, const int* v,
             const TReal* ay, int lday);

    template<
      template<typename R>
      class TPreorder
      >
    DcaState(int size, int mark, const int* v,
             const TReal* rz, int ldrz,
             const TPreorder<TReal>& preo);

    template<
      template<typename R>
      class TPreorder
      >
    DcaState(int m, int size, int mark, const int* v,
             const TReal* ay, int lday,
             const TPreorder<TReal>& preo);

    void
    dropColumn(int mark);

    void
    nextNode();

    template<
      template<typename R>
      class TPreorder
      >
    void
    nextNode(const TPreorder<TReal>& preo);

    bool
    isFinal() const;

    int
    currentSize() const;

    int
    currentMark() const;

    TReal
    minBound() const;

    TReal
    minBound(int size, const TReal* tau) const;

    template<
      template<typename R>
      class TCriterion
    >
    TReal
    minBound(int mark, const TCriterion<TReal>& crit) const;

    template<
      template<typename R>
      class TCriterion
    >
    void
    reportSubleading(DcaTable<TReal,TCriterion>& table) const;

    int
    nodeCount() const;

  };


}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/dca_state.cc"
#endif
