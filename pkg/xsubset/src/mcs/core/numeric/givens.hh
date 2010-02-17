/**
 * @file givens.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_NUMERIC_GIVENS_HH
#define MCS_CORE_NUMERIC_GIVENS_HH


#define RANGE range<Size>
#define VECTOR vector<value_type, Alloc>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Size>
      class range;


      template<typename Value,
	       typename Alloc>
      class vector;


      template<typename Value>
      class givens
      {

      public:

	typedef Value value_type;

	typedef value_type& reference;

	typedef const value_type& const_reference;

	typedef value_type* pointer;

	typedef const value_type* const_pointer;

	typedef size_t size_type;


      private:

	value_type r_;

	value_type c_;

	value_type s_;


      public:

	givens();

	givens(value_type dx, value_type dy);


	value_type
	r() const;

	value_type
	c() const;

	value_type
	s() const;

	void
	gen(value_type dx, value_type dy);

	void
	rot(reference x, reference y);

	void
	rot(size_type n,
	    reference x, size_type incx,
	    reference y, size_type incy);

	void
	rot(size_type n,
	    const_reference x, size_type incx,
	    const_reference y, size_type incy,
	    reference xx, size_type incxx,
	    reference yy, size_type incyy);

	template<typename Alloc>
	void
	rot(VECTOR&& x, VECTOR&& y);

	template<typename Alloc>
	void
	rot(const VECTOR& x, const VECTOR& y,
	    VECTOR&& xx, VECTOR&& yy);


	template<typename Alloc>
	static
	void
	zero(VECTOR&& x, VECTOR&& y);

	template<typename Size,
		 typename Alloc>
	static
	void
	zero(RANGE r, VECTOR&& x, VECTOR&& y);

	template<typename Alloc>
	static
	void
	zero(const VECTOR& x, const VECTOR& y,
	     VECTOR&& xx, VECTOR&& yy);

	template<typename Size,
		 typename Alloc>
	static
	void
	zero(RANGE r, const VECTOR& x, const VECTOR& y,
	     RANGE rr, VECTOR&& xx, VECTOR&& yy);
		 
      };


    }

  }

}


#undef RANGE
#undef VECTOR


#include "givens.cc"
#endif
