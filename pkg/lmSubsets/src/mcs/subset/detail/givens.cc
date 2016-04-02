#ifndef MCS_SUBSET_DETAIL_GIVENS_CC
#define MCS_SUBSET_DETAIL_GIVENS_CC



#include "mcs/subset/detail/givens.hh"

#include <cmath>  // std::abs, std::sqrt, std::copysign



namespace mcs    {
namespace subset {
namespace detail {



  template<typename TReal>
  Givens<TReal>::Givens() :
    r_{0}, c_{0}, s_{0}
  {
  }


  template<typename TReal>
  Givens<TReal>::Givens(const TReal dx, const TReal dy)
  {
    gen(dx, dy);
  }


  template<typename TReal>
  TReal
  Givens<TReal>::r() const
  {
    return r_;
  }


  template<typename TReal>
  TReal
  Givens<TReal>::c() const
  {
    return c_;
  }


  template<typename TReal>
  TReal
  Givens<TReal>::s() const
  {
    return s_;
  }


  template<typename TReal>
  void
  Givens<TReal>::gen(const TReal dx, const TReal dy)
  {
    // Source: Wikipedia
    if (dy == 0)
      {
        c_ = std::copysign(1, dx);
        s_ = 0;
        r_ = std::abs(dx);
      }
    else if (dx == 0)
      {
        c_ = 0;
        s_ = std::copysign(1, dy);
        r_ = std::abs(dy);
      }
    else if (std::abs(dy) > std::abs(dx))
      {
        TReal t = dx / dy;
        TReal u = std::copysign(std::sqrt(1 + t * t), dy);

        s_ = 1 / u;
        c_ = s_ * t;
        r_ = dy * u;
      }
    else
      {
        TReal t = dy / dx;
        TReal u = std::copysign(std::sqrt(1 + t * t), dx);

        c_ = 1 / u;
        s_ = c_ * t;
        r_ = dx * u;
      }
  }


  template<typename TReal>
  void
  Givens<TReal>::rot(const int n, TReal* x, const int incx,
                                  TReal* y, const int incy)
  {
    for (int i = 0; i < n; ++i)
      {
        TReal  t =  c_ * (*x) + s_ * (*y);
              *y = -s_ * (*x) + c_ * (*y);
              *x =  t;

        x += incx;
        y += incy;
      }
  }


  template<typename TReal>
  void
  Givens<TReal>::rot(const int n, const TReal* x , const int incx ,
                                  const TReal* y , const int incy ,
                                        TReal* xx, const int incxx,
                                        TReal* yy, const int incyy)
  {
    for (int i = 0; i < n; ++i)
      {
        TReal   t =  c_ * (*x) + s_ * (*y);
              *yy = -s_ * (*x) + c_ * (*y);
              *xx =  t;

        x  += incx ;  y  += incy ;
        xx += incxx;  yy += incyy;
      }
  }


  template<typename TReal>
  void
  Givens<TReal>::zero(const int n, TReal* x, const int incx,
                                   TReal* y, const int incy)
  {
    if (n > 0)
      {
	Givens g(*x, *y);
	g.rot(n - 1, x + incx, incx, y + incx, incy);

	*x = g.r_;
	*y = 0;
      }
  }


  template<typename TReal>
  void
  Givens<TReal>::zero(const int n, const TReal* x , const int incx ,
                                   const TReal* y , const int incy ,
                                         TReal* xx, const int incxx,
                                         TReal* yy, const int incyy)
  {
    if (n > 0)
      {
	Givens g(*x, *y);
	g.rot(n - 1,  x + incx , incx ,  y + incy , incy ,
                     xx + incxx, incxx, yy + incyy, incyy);

	*xx = g.r_;
	*yy = 0;
      }
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif
