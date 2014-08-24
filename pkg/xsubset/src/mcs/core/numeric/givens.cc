/**
 * @file givens.cc
 */
#ifndef MCS_CORE_NUMERIC_GIVENS_CC
#define MCS_CORE_NUMERIC_GIVENS_CC


#include <cmath>

#include "../../mcs.hh"

#include "slice.hh"
#include "vector.hh"
#include "givens.hh"


#define GIVENS givens<Value>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value>
      GIVENS::givens() :
        r_(0), s_(0), c_(0)
      {
      }


      template<typename Value>
      GIVENS::givens(const Value dx,
		     const Value dy)
      {
	gen(dx, dy);
      }


      template<typename Value>
      Value
      GIVENS::r() const
      {
	return r_;
      }


      template<typename Value>
      Value
      GIVENS::c() const
      {
	return c_;
      }


      template<typename Value>
      Value
      GIVENS::s() const
      {
	return s_;
      }


      template<typename Value>
      void
      GIVENS::gen(const Value dx, const Value dy)
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
	    Value t = dx / dy;
	    Value u = std::copysign(std::sqrt(1 + t * t), dy);
	    s_ = 1 / u;
	    c_ = s_ * t;
	    r_ = dy * u;
	  }
	else
	  {
	    Value t = dy / dx;
	    Value u = std::copysign(std::sqrt(1 + t * t), dx);
	    c_ = 1 / u;
	    s_ = c_ * t;
	    r_ = dx * u;
	  }
      }


      template<typename Value>
      void
      GIVENS::rot(Value& x, Value& y)
      {
        Value& t =  c_ * x + s_ * y;
	y        = -s_ * x + c_ * y;
	x        =  t;
      }


      template<typename Value>
      template<typename Size>
      void
      GIVENS::rot(const Size n,
		  Value& x, const Size incx,
		  Value& y, const Size incy)
      {
        MCS_ASSERT(n >= 0, "invalid argument (givens::rot)");
        MCS_ASSERT(incx > 0, "invalid argument (givens::rot)");
        MCS_ASSERT(incy > 0, "invalid argument (givens::rot)");

	Value* posx = &x;
	Value* posy = &y;

        for (Size i = 0; i < n; ++i)
	  {
	    Value t =  c_ * (*posx) + s_ * (*posy);
	    *posy   = -s_ * (*posx) + c_ * (*posy);
	    *posx   =  t;

            posx += incx;
            posy += incy;
	  }
      }


      template<typename Value>
      template<typename Size>
      void
      GIVENS::rot(const Size n,
		  const Value& x, const Size incx,
		  const Value& y, const Size incy,
		  Value& xx, const Size incxx,
		  Value& yy, const Size incyy)
      {
        MCS_ASSERT(n >= 0, "invalid argument (givens::rot)");
        MCS_ASSERT(incx > 0, "invalid argument (givens::rot)");
        MCS_ASSERT(incy > 0, "invalid argument (givens::rot)");
        MCS_ASSERT(incxx > 0, "invalid argument (givens::rot)");
        MCS_ASSERT(incyy > 0, "invalid argument (givens::rot)");

        const Value* posx = &x;
	const Value* posy = &y;
	Value* posxx = &xx;
	Value* posyy = &yy;

        for (Size i = 0; i < n; ++i)
	  {
	    Value t =  c_ * (*posx) + s_ * (*posy);
	    *posyy  = -s_ * (*posx) + c_ * (*posy);
	    *posxx  = t;

            posx += incx;    posy += incy;
            posxx += incxx;  posyy += incyy;
	  }
      }


      template<typename Value>
      template<typename Size,
               template<typename V,
                        typename S>
               class Derived1,
               template<typename V,
                        typename S>
               class Derived2>
      void
      GIVENS::rot(vector_base<Value, Size, Derived1>&& x,
                  vector_base<Value, Size, Derived2>&& y)
      {
        MCS_ASSERT(x.len() == y.len(), "invalid argument (givens::rot)");

	rot(x.len(), x(0), x.inc(), y(0), y.inc());
      }


      template<typename Value>
      template<typename Size,
               template<typename V,
                        typename S>
               class Derived1,
               template<typename V,
                        typename S>
               class Derived2,
               template<typename V,
                        typename S>
               class Derived3,
               template<typename V,
                        typename S>
               class Derived4>
      void
      GIVENS::rot(const vector_base<Value, Size, Derived1>& x,
                  const vector_base<Value, Size, Derived2>& y,
                  vector_base<Value, Size, Derived3>&& xx,
                  vector_base<Value, Size, Derived4>&& yy)
      {
        MCS_ASSERT(x.len() == y.len(), "invalid argument (givens::rot)");
        MCS_ASSERT(xx.len() == x.len(), "invalid argument (givens::rot)");
        MCS_ASSERT(yy.len() == y.len(), "invalid argument (givens::rot)");

	rot(x.len(), x(0), x.inc(), y(0), y.inc(),
	    xx(0), xx.inc(), yy(0), yy.inc());
      }


      template<typename Value>
      template<typename Size,
               template<typename V,
                        typename S>
               class Derived1,
               template<typename V,
                        typename S>
               class Derived2>
      void
      GIVENS::zero(vector_base<Value, Size, Derived1>&& x,
                   vector_base<Value, Size, Derived2>&& y)
      {
        //std::cout << "givens::zero(vector_base&&, vector_base&&)" << std::endl;

        MCS_ASSERT(x.len() == y.len(), "invalid argument (givens::zero)");

        const Size n = x.len();

        if (n == 0)
          return;

	givens<Value> g(x(0), y(0));
        if (n > 1)
          g.rot(n - 1, x(1), x.inc(), y(1), y.inc());
	x(0) = g.r_;
	y(0) = 0;
      }


      template<typename Value>
      template<typename Size,
               template<typename V,
                        typename S>
               class Derived1,
               template<typename V,
                        typename S>
               class Derived2>
      void
      GIVENS::zero(const slice<Size>& s,
                   vector_base<Value, Size, Derived1>&& x,
                   vector_base<Value, Size, Derived2>&& y)
      {
        //std::cout << "givens::zero(slice, vector_base&&, vector_base&&)" << std::endl;

        MCS_ASSERT(s.pos() + s.len() <= x.len(),
                   "invalid argument (givens::zero)");
        MCS_ASSERT(s.pos() + s.len() <= y.len(),
                   "invalid argument (givens::zero)");

	zero(x(s), y(s));
      }
      

      template<typename Value>
      template<typename Size,
               template<typename V,
                        typename S>
               class Derived1,
               template<typename V,
                        typename S>
               class Derived2,
               template<typename V,
                        typename S>
               class Derived3,
               template<typename V,
                        typename S>
               class Derived4>
      void
      GIVENS::zero(const vector_base<Value, Size, Derived1>& x,
                   const vector_base<Value, Size, Derived2>& y,
                   vector_base<Value, Size, Derived3>&& xx,
                   vector_base<Value, Size, Derived4>&& yy)
      {
        MCS_ASSERT(x.len() == y.len(), "invalid argument (givens::rot)");
        MCS_ASSERT(xx.len() == x.len(), "invalid argument (givens::rot)");
        MCS_ASSERT(yy.len() == y.len(), "invalid argument (givens::rot)");

        const Size n = x.len();

        if (n == 0)
          return;

	givens<Value> g(x(0), y(0));
        if (n > 1)
          g.rot(n - 1, x(1), x.inc(), y(1), y.inc(),
                xx(1), xx.inc(), yy(1), yy.inc());
	xx(0) = g.r_;
	yy(0) = 0;
      }
    

      template<typename Value>
      template<typename Size,
               template<typename V,
                        typename S>
               class Derived1,
               template<typename V,
                        typename S>
               class Derived2,
               template<typename V,
                        typename S>
               class Derived3,
               template<typename V,
                        typename S>
               class Derived4>
      void
      GIVENS::zero(const slice<Size>& s,
                   const vector_base<Value, Size, Derived1>& x,
                   const vector_base<Value, Size, Derived2>& y,
                   const slice<Size>& ss,
                   vector_base<Value, Size, Derived3>&& xx,
                   vector_base<Value, Size, Derived4>&& yy)
      {
        MCS_ASSERT(s.pos() + s.len() <= x.len(),
                   "invalid argument (givens::zero)");
        MCS_ASSERT(s.pos() + s.len() <= y.len(),
                   "invalid argument (givens::zero)");
        MCS_ASSERT(ss.pos() + ss.len() <= xx.len(),
                   "invalid argument (givens::zero)");
        MCS_ASSERT(ss.pos() + ss.len() <= yy.len(),
                   "invalid argument (givens::zero)");

	zero(x(s), y(s), xx(ss), yy(ss));
      }


    }

  }

}


#endif
