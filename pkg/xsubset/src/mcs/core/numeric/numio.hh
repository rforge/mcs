/**
 * @file numio.hh
 */
#ifndef MCS_CORE_NUMERIC_NUMIO_HH
#define MCS_CORE_NUMERIC_NUMIO_HH


#include <iostream>

#include "vector.hh"
#include "matrix.hh"


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
                  vector_base<Value, Size, Derived>&& x);

   
      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      std::ostream&
      operator <<(std::ostream& os,
                  const vector_base<Value, Size, Derived>& x);


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      std::istream&
      operator >>(std::istream& is,
                  matrix_base<Value, Size, Derived>&& a);


      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      std::ostream&
      operator <<(std::ostream& os,
                  const matrix_base<Value, Size, Derived>& a);


    }

  }

}


#include "numio.cc"
#endif
