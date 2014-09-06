#ifndef MCS_SUBSET_DCA_HH
#define MCS_SUBSET_DCA_HH


#include "mcs/subset/detail/math.hh"
#include "mcs/subset/detail/givens.hh"
#include "mcs/subset/detail/lapack.hh"
#include "mcs/subset/detail/qrz.hh"
#include "mcs/subset/detail/dca_table.hh"
#include "mcs/subset/detail/dca_state.hh"

#include "mcs/subset/detail/dca.hh"



namespace mcs    {
namespace subset {



  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  int dca(const int m, const int size, const int mark, const int nbest,
          const int* const v, const TReal* const ay, const int lday,
	  int* const sIndex, TReal* const sRss, TReal* const sCrit,
          int* const sSize, int* const s, const TCriterion<TReal>& c)
  {
    using namespace detail;

    DcaState<TReal>            state(m, size, mark, v, ay, lday);
    DcaTable<TReal,TCriterion> table(size, nbest, sIndex, sRss, sCrit, sSize, s, c);

    return dca(state, table);
  }


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  int dca(const int size, const int mark, const int nbest,
          const int* const v, const TReal* const rz, const int ldrz,
	  int* const sIndex, TReal* const sRss, TReal* const sCrit,
          int* const sSize, int* const s, const TCriterion<TReal>& c)
  {
    using namespace detail;

    DcaState<TReal>            state(size, mark, v, rz, ldrz);
    DcaTable<TReal,TCriterion> table(size, nbest, sIndex, sRss, sCrit, sSize, s, c);

    return dca(state, table);
  }


  template<typename TReal>
  int dca(const int m, const int size, const int mark, const int nbest,
          const int* const v, const TReal* const ay, const int lday,
	  int* const sIndex, TReal* const sRss, int* const s)
  {
    using namespace detail;

    DcaState<TReal>                state(m, size, mark, v, ay, lday);
    DcaTable<TReal,Criteria::None> table(size, mark, nbest, sIndex, sRss, s);

    return dca(state, table);
  }


  template<typename TReal>
  int dca(const int size, const int mark, const int nbest,
          const int* const v, const TReal* const rz, const int ldrz,
	  int* const sIndex, TReal* const sRss, int* const s)
  {
    using namespace detail;

    DcaState<TReal>                state(size, mark, v, rz, ldrz);
    DcaTable<TReal,Criteria::None> table(size, mark, nbest, sIndex, sRss, s);

    return dca(state, table);
  }



}  // namespace subset
}  // namespace mcs



#endif
