/**
 * @file mcs/mcs.hh
 *
 * @brief General definitions for the MCS namespace.
 *
 * @author Marc Hofmann
 */


#ifndef _MCS_HH_
#define _MCS_HH_


#include <cassert>


#ifndef MCS_INLINE
#define MCS_INLINE inline
#endif


#ifndef MCS_CHECK_INTERRUPT
#define MCS_CHECK_INTERRUPT
#endif


//#define NDEBUG  // disable asserts

#define MCS_ASSERT(expr) \
  assert(expr)



/**
 * @brief MCS namespace.
 */
namespace mcs
{
}


#endif
