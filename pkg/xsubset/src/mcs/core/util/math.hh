/**
 * @file math.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_UTIL_MATH_HH
#define MCS_CORE_UTIL_MATH_HH


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
	sqr(Value val)
	{
	  return val * val;
	}


	template<typename Value>
	struct constants
	{

	  typedef Value value_type;

	  static const value_type pi = PI;

	};


      }

    }

  }

}


#undef PI


#endif
