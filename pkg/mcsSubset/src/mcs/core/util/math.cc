/**
 * @file math.cc
 */
#ifndef MCS_CORE_UTIL_MATH_CC
#define MCS_CORE_UTIL_MATH_CC


#include "math.hh"


#define PI 3.14159265358979323846264338327


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
	sign(const Value dx)
	{
          return (dx > 0) - (dx < 0);
	}


	template<typename Value>
	Value
	sqr(Value val)
	{
	  return val * val;
	}


	template<typename Value>
	const Value constants<Value>::pi = PI;


      }

    }

  }

}


#undef PI


#endif
