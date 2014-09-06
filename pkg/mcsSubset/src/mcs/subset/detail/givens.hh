#ifndef MCS_CORE_SUBSET_GIVENS_HH
#define MCS_CORE_SUBSET_GIVENS_HH



namespace mcs    {
namespace subset {
namespace detail {



  template<typename TReal>
  class Givens
  {

  private:

    TReal r_;
    TReal c_;
    TReal s_;


  public:

    Givens();
    Givens(TReal dx, TReal dy);


    TReal r() const;
    TReal c() const;
    TReal s() const;


    void
    gen(TReal dx, TReal dy);


    void
    rot(int n, TReal* x, int incx,
               TReal* y, int incy);


    void
    rot(int n, const TReal* x , int incx , const TReal* y , int incy ,
                     TReal* xx, int incxx,       TReal* yy, int incyy);


    static
    void
    zero(int n, TReal* x, int incx,
                TReal* y, int incy);


    static
    void
    zero(int n, const TReal* x , int incx , const TReal* y , int incy ,
                      TReal* xx, int incxx,       TReal* yy, int incyy);

  };



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/givens.cc"
#endif
