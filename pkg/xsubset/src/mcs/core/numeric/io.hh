/**
 * @file io.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_NUMERIC_IO_HH
#define MCS_CORE_NUMERIC_IO_HH


#include <iostream>


#define VECTOR vector<Value, Alloc>
#define MATRIX matrix<Value, Alloc>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
	       typename Alloc>
      class vector;


      template<typename Value,
	       typename Alloc>
      class matrix;


      template<typename Value,
	       typename Alloc>
      std::istream&
      operator >>(std::istream& is, VECTOR&& vec);

   
      template<typename Value,
	       typename Alloc>
      std::ostream&
      operator <<(std::ostream& os, const VECTOR& vec);


      template<typename Value,
	       typename Alloc>
      std::istream&
      operator >>(std::istream& is, MATRIX&& mat);


      template<typename Value,
	       typename Alloc>
      std::ostream&
      operator <<(std::ostream& os, const MATRIX& mat);


    }

  }

}


#undef VECTOR
#undef MATRIX


#include "io.cc"
#endif
