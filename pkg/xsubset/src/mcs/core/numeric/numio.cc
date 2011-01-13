/**
 * @file numio.cc
 */
#ifndef MCS_CORE_NUMERIC_NUMIO_CC
#define MCS_CORE_NUMERIC_NUMIO_CC


#include <iostream>

#include "vector.hh"
#include "matrix.hh"
#include "numio.hh"


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      std::istream&
      operator >>(std::istream& is,
                  vector_base<Value, Size, Derived>&& x)
      {
	for (Size i = 0; i < x.len(); ++i)
	  {
	    is >> x(i);
	  }

	return is;
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      std::ostream&
      operator <<(std::ostream& os,
                  const vector_base<Value, Size, Derived>& x)
      {
	for (Size i = 0; i < x.len(); ++i)
	  {
	    os << " " << x(i);
	  }

	return os;
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      std::istream&
      operator >>(std::istream& is,
                  matrix_base<Value, Size, Derived>&& a)
      {
	for (Size i = 0; i < a.nrow(); ++i)
	  {
	    for (Size j = 0; j < a.ncol(); ++j)
	      {
		is >> a(i, j);
	      }
	  }

	return is;
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      std::ostream&
      operator <<(std::ostream& os,
                  const matrix_base<Value, Size, Derived>& a)
      {
	for (Size i = 0; i < a.nrow(); ++i)
	  {
	    for (Size j = 0; j < a.ncol(); ++j)
	      {
		os << " " << a(i, j);
	      }

	    os << std::endl;
	  }

	return os;
      }


    }

  }

}


#endif
