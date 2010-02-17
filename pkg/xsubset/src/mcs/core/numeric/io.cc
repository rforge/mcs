/**
 * @file io.cc
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_NUMERIC_IO_CC
#define MCS_CORE_NUMERIC_IO_CC


#include <iostream>

#include "vector.hh"
#include "matrix.hh"
#include "io.hh"


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
      std::istream&
      operator >>(std::istream& is, VECTOR&& vec)
      {
	for (int i = 0; i < vec.len(); ++i)
	  {
	    is >> vec(i);
	  }

	return is;
      }


      template<typename Value,
	       typename Alloc>
      std::ostream&
      operator <<(std::ostream& os, const VECTOR& vec)
      {
	for (int i = 0; i < vec.len(); ++i)
	  {
	    os << " " << vec(i);
	  }
	os << std::endl;

	return os;
      }


      template<typename Value,
	       typename Alloc>
      std::istream&
      operator >>(std::istream& is, MATRIX&& mat)
      {
	for (int i = 0; i < mat.nrow(); ++i)
	  {
	    for (int j = 0; j < mat.ncol(); ++j)
	      {
		is >> mat(i, j);
	      }
	  }

	return is;
      }


      template<typename Value,
	       typename Alloc>
      std::ostream&
      operator <<(std::ostream& os, const MATRIX& mat)
      {
	for (int i = 0; i < mat.nrow(); ++i)
	  {
	    for (int j = 0; j < mat.ncol(); ++j)
	      {
		os << " " << mat(i, j);
	      }

	    os << std::endl;
	  }

	return os;
      }


    }

  }

}


#undef VECTOR
#undef MATRIX


#endif
