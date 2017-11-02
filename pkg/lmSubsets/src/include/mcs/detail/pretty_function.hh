#ifndef MCS_DETAIL_PRETTY_FUNCTION_HH
#define MCS_DETAIL_PRETTY_FUNCTION_HH

// see also:  http://www.boost.org/libs/utility/current_function.html



namespace mcs    {
namespace detail {



inline void
aux_pretty_function_()
{
#if   defined( __GNUC__)
#  define MCS_PRETTY_FUNCTION __PRETTY_FUNCTION__
#elif defined(__FUNCSIG__)
#  define MCS_PRETTY_FUNCTION __FUNCSIG__
#else
#  define MCS_PRETTY_FUNCTION __func__
#endif
}



}  // end namespace mcs
}  // end namespace detail



#endif
