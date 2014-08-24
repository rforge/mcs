/**
 * @file givens.hh
 */
#ifndef MCS_CORE_NUMERIC_GIVENS_HH
#define MCS_CORE_NUMERIC_GIVENS_HH


#include "slice.hh"
#include "vector.hh"
#include "matrix.hh"


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value>
      class givens
      {

      private:

	Value r_;

	Value c_;

	Value s_;


      public:

	givens();

	givens(Value dx, Value dy);


        Value
	r() const;

        Value
	c() const;

        Value
	s() const;

	void
	gen(Value dx, Value dy);

	void
	rot(Value& x, Value& y);

        template<typename Size>
	void
	rot(Size n,
	    Value& x, Size incx,
	    Value& y, Size incy);

        template<typename Size>
	void
	rot(Size n,
	    const Value& x, Size incx,
	    const Value& y, Size incy,
	    Value& xx, Size incxx,
	    Value& yy, Size incyy);
        
	template<typename Size,
                 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2>
	void
        rot(vector_base<Value, Size, Derived1>&& x,
            vector_base<Value, Size, Derived2>&& y);

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
	rot(const vector_base<Value, Size, Derived1>& x,
            const vector_base<Value, Size, Derived2>& y,
	    vector_base<Value, Size, Derived3>&& xx,
            vector_base<Value, Size, Derived4>&& yy);
        
        template<typename Size,
                 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2>
	static
	void
	zero(vector_base<Value, Size, Derived1>&& x,
             vector_base<Value, Size, Derived2>&& y);

        template<typename Size,
                 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2>
	static
	void
	zero(const slice<Size>& s,
             vector_base<Value, Size, Derived1>&& x,
             vector_base<Value, Size, Derived2>&& y);

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
	static
	void
	zero(const vector_base<Value, Size, Derived1>& x,
             const vector_base<Value, Size, Derived2>& y,
	     vector_base<Value, Size, Derived3>&& xx,
             vector_base<Value, Size, Derived4>&& yy);
        
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
	static
	void
	zero(const slice<Size>& s,
             const vector_base<Value, Size, Derived1>& x,
             const vector_base<Value, Size, Derived2>& y,
             const slice<Size>& ss,
	     vector_base<Value, Size, Derived3>&& xx,
             vector_base<Value, Size, Derived4>&& yy);
		 
      };


    }

  }

}


#include "givens.cc"
#endif
