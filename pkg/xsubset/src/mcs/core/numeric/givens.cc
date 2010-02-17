/**
 * @file givens.cc
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_NUMERIC_GIVENS_CC
#define MCS_CORE_NUMERIC_GIVENS_CC


#include <cmath>

#include "vector.hh"
#include "givens.hh"


#define RANGE range<Size>
#define VECTOR vector<value_type, Alloc>
#define GIVENS givens<Value>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {

      namespace detail
      {


	float
	sign(float dx)
	{
	  if (dx > 0.0)
	    {
	      return 1.0;
	    }
	  else if (dx < 0.0)
	    {
	      return -1.0;
	    }
	  else
	    {
	      return 0.0;
	    }
	}


	double
	sign(double dx)
	{
	  if (dx > 0.0L)
	    {
	      return 1.0L;
	    }
	  else if (dx < 0.0L)
	    {
	      return -1.0L;
	    }
	  else
	    {
	      return 0.0L;
	    }
	}


	float
	copysign(float dx, float dy)
	{
	  return (dx * sign(dy));
	}


	double
	copysign(double dx, double dy)
	{
	  return (dx * sign(dy));
	}


      }


      template<typename Value>
      GIVENS::givens() : r_(0), s_(0), c_(0)
      {
      }


      template<typename Value>
      GIVENS::givens(const value_type dx,
		     const value_type dy)
      {
	gen(dx, dy);
      }


      template<typename Value>
      typename GIVENS::value_type
      GIVENS::r() const
      {
	return r_;
      }


      template<typename Value>
      typename GIVENS::value_type
      GIVENS::c() const
      {
	return c_;
      }


      template<typename Value>
      typename GIVENS::value_type
      GIVENS::s() const
      {
	return s_;
      }


      template<typename Value>
      void
      GIVENS::gen(const value_type dx, const value_type dy)
      {
	if (dy == 0.0L)
	  {
	    c_ = detail::copysign(1.0L, dx);
	    s_ = 0.0L;
	    r_ = std::abs(dx);
	  }
	else if (dx == 0)
	  {
	    c_ = 0.0L;
	    s_ = detail::copysign(1.0L, dy);
	    r_ = std::abs(dy);
	  }
	else if (std::abs(dy) > std::abs(dx))
	  {
	    value_type t = dx / dy;
	    value_type u = detail::copysign(std::sqrt(1 + t * t), dy);
	    s_ = 1 / u;
	    c_ = s_ * t;
	    r_ = dy * u;
	  }
	else
	  {
	    value_type t = dy / dx;
	    value_type u = detail::copysign(std::sqrt(1 + t * t), dx);
	    c_ = 1 / u;
	    s_ = c_ * t;
	    r_ = dx * u;
	  }
      }


      template<typename Value>
      void
      GIVENS::rot(reference x, reference y)
      {
	value_type t =  c_ * x + s_ * y;
	y            = -s_ * x + c_ * y;
	x            =  t;
      }


      template<typename Value>
      void
      GIVENS::rot(const size_type n,
		  reference x, const size_type incx,
		  reference y, const size_type incy)
      {
	pointer posx = &x;
	pointer posy = &y;
	const_pointer endx = posx + n * incx;

	for (; posx < endx; posx += incx, posy += incy)
	  {
	    value_type t =  c_ * (*posx) + s_ * (*posy);
	    *posy        = -s_ * (*posx) + c_ * (*posy);
	    *posx        =  t;
	  }
      }


      template<typename Value>
      void
      GIVENS::rot(const size_type n,
		  const_reference x, const size_type incx,
		  const_reference y, const size_type incy,
		  reference xx, const size_type incxx,
		  reference yy, const size_type incyy)
      {
	const_pointer posx = &x;
	const_pointer posy = &y;
	pointer posxx = &xx;
	pointer posyy = &yy;
	const_pointer const endx = posx + n * incx;

	for (;
	     posx < endx;
	     posx += incx, posy += incy, posxx += incxx, posyy += incyy)
	  {
	    value_type t =  c_ * (*posx) + s_ * (*posy);
	    *posyy       = -s_ * (*posx) + c_ * (*posy);
	    *posxx       = t;
	  }
      }


      template<typename Value>
      template<typename Alloc>
      void
      GIVENS::rot(VECTOR&& x, VECTOR&& y)
      {
	rot(x.len(), x(0), x.inc(), y(0), y.inc());
      }


      template<typename Value>
      template<typename Alloc>
      void
      GIVENS::rot(const VECTOR& x, const VECTOR& y,
		  VECTOR&& xx, VECTOR&& yy)
      {
	rot(x.len(), x(0), x.inc(), y(0), y.inc(),
	    xx(0), xx.inc(), yy(0), yy.inc());
      }


      template<typename Value>
      template<typename Alloc>
      void
      GIVENS::zero(VECTOR&& x, VECTOR&& y)
      {
	givens g(x(0), y(0));
	g.rot(x.len() - 1, x(1), x.inc(), y(1), y.inc());
	x(0) = g.r_;
	y(0) = 0;
      }


      template<typename Value>
      template<typename Size,
	       typename Alloc>
      void
      GIVENS::zero(RANGE r, VECTOR&& x, VECTOR&& y)
      {
	zero(x(r), y(r));
      }


      template<typename Value>
      template<typename Alloc>
      void
      GIVENS::zero(const VECTOR& x, const VECTOR& y,
		   VECTOR&& xx, VECTOR&& yy)
      {
	givens g(x(0), y(0));
	g.rot(x.len() - 1, x(1), x.inc(), y(1), y.inc(),
	      xx(1), xx.inc(), yy(1), yy.inc());
	xx(0) = g.r_;
	yy(0) = 0;
      }


      template<typename Value>
      template<typename Size,
	       typename Alloc>
      void
      GIVENS::zero(RANGE r, const VECTOR& x, const VECTOR& y,
		   RANGE rr, VECTOR&& xx, VECTOR&& yy)
      {
	zero(x(r), y(r), xx(rr), yy(rr));
      }


    }

  }

}


#undef RANGE
#undef VECTOR
#undef GIVENS


#endif
