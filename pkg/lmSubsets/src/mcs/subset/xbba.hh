#ifndef MCS_SUBSET_XBBA_HH
#define MCS_SUBSET_XBBA_HH


#include "mcs/subset/detail/algo.hh"
#include "mcs/subset/detail/math.hh"
#include "mcs/subset/detail/givens.hh"
#include "mcs/subset/detail/lapack.hh"
#include "mcs/subset/detail/qrz.hh"
#include "mcs/subset/detail/preorder.hh"
#include "mcs/subset/detail/dca_table.hh"
#include "mcs/subset/detail/dca_state.hh"

#include "mcs/subset/detail/pbba.hh"

#include <cstring>



namespace mcs    {
namespace subset {



  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  int
  xbba(const char* const algo,
       const int m, const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const ay, const int lday,
       int* const sIndex, TReal* const sRss, TReal* const sCrit,
       int* const sSize, int* const s, const TCriterion<TReal>& c,
       int* const info)
  {
    using namespace detail;


    Preorder::Full<TReal> p1(size, pmin);

    DcaState<TReal>            state(m, size, mark, v, ay, lday, p1);
    DcaTable<TReal,TCriterion> table(size, nbest, sIndex, sRss, sCrit, sSize, s, c);


    if (std::strcmp(algo, "xbba1")) {
      Preorder::Single1<TReal> p(size, pmin);

      return pbba(state, table, c, p);
    } else if (std::strcmp(algo, "xbba2")) {
      Preorder::Single2<TReal> p(size, pmin);

      return pbba(state, table, c, p);
    } else {
      *info = 1;

      return -1;
    }
  }


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  int
  xbba(const char* const algo,
       const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const rz, const int ldrz,
       int* const sIndex, TReal* const sRss, TReal* const sCrit,
       int* const sSize, int* const s, const TCriterion<TReal>& c,
       int* const info)
  {
    using namespace detail;


    Preorder::Full<TReal> p1(size, pmin);

    DcaState<TReal>            state(size, mark, v, rz, ldrz, p1);
    DcaTable<TReal,TCriterion> table(size, nbest, sIndex, sRss, sCrit, sSize, s, c);


    if (std::strcmp(algo, "xbba1")) {
      Preorder::Single1<TReal> p (size, pmin);

      return pbba(state, table, c, p);
    } else if (std::strcmp(algo, "xbba2")) {
      Preorder::Single2<TReal> p (size, pmin);

      return pbba(state, table, c, p);
    } else {
      *info = 1;

      return -1;
    }
  }


  template<typename TReal>
  int
  xbba(const char* const algo,
       const int m, const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const ay, const int lday,
       int* const sIndex, TReal* const sRss, int* const s,
       int* const info)
  {
    using namespace detail;


    Preorder::Full<TReal> p1(size, pmin);

    DcaState<TReal>                 state(m, size, mark, v, ay, lday, p1);
    DcaTable<TReal, Criteria::None> table(size, mark, nbest, sIndex, sRss, s);


    if (std::strcmp(algo, "xbba1")) {
      Preorder::Single1<TReal> p (size, pmin);

      return pbba(state, table, p);
    } else if (std::strcmp(algo, "xbba2")) {
      Preorder::Single2<TReal> p (size, pmin);

      return pbba(state, table, p);
    } else {
      *info = 1;

      return -1;
    }
  }


  template<typename TReal>
  int
  xbba(const char* const algo,
       const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const rz, const int ldrz,
       int* const sIndex, TReal* const sRss, int* const s,
       int* const info)
  {
    using namespace detail;


    Preorder::Full<TReal> p1(size, pmin);

    DcaState<TReal>                 state(size, mark, v, rz, ldrz, p1);
    DcaTable<TReal, Criteria::None> table(size, mark, nbest, sIndex, sRss, s);


    if (std::strcmp(algo, "xbba1")) {
      Preorder::Single1<TReal> p (size, pmin);

      return pbba(state, table, p);
    } else if (std::strcmp(algo, "xbba2")) {
      Preorder::Single2<TReal> p (size, pmin);

      return pbba(state, table, p);
    } else {
      *info = 1;

      return -1;
    }
  }



}  // namespace subset
}  // namespace mcs



#endif
