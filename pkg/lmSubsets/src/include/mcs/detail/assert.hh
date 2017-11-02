#ifndef MCS_DETAIL_ASSERT_HH
#define MCS_DETAIL_ASSERT_HH

// see also:  http://www.boost.org/libs/assert/assert.html



#include <cstdlib>  // std::abort
#include <iostream>  // std::cerr, std::endl



#include "mcs/detail/pretty_function.hh"



#define MCS_ASSERT(expr, msg)                                           \
    ((expr)?                                                            \
     ((void) 0) :                                                       \
     mcs::detail::assertion_failed(#expr, #msg, MCS_PRETTY_FUNCTION,    \
                                   __FILE__, __LINE__))



namespace mcs    {
namespace detail {



void
assertion_failed(const char* expr, const char* msg, const char* func,
                 const char* file, const long line)
{
  std::cerr << file << ":" << line << ":  " << "Assertion failed (in '" << func
            << "'):  " << msg << " (" << expr << ")" << std::endl;
  std::abort();
}



}  // end namespace mcs
}  // end namespace detail



#endif
