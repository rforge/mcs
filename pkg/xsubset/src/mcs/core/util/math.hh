/**
 * @file math.hh
 */
#ifndef MCS_CORE_UTIL_MATH_HH
#define MCS_CORE_UTIL_MATH_HH


namespace mcs
{

  namespace core
  {

    namespace util
    {

      namespace math
      {


        template<typename Value>
        Value
	sign(const Value dx);


	template<typename Value>
	Value
	sqr(Value val);


	template<typename Value>
	struct constants
	{

	  static const Value pi;

	};


      }

    }

  }

}


#include "math.cc"
#endif
