/**
 * @file mcs.hh
 *
 * @brief General definitions for the MCS namespace.
 */
#ifndef MCS_HH
#define MCS_HH


#include <cassert>
#include <iostream>


#ifndef MCS_INLINE
#define MCS_INLINE inline
#endif


#ifndef MCS_CHECK_INTERRUPT
#define MCS_CHECK_INTERRUPT
#endif


/**
 * Disable degug.
 */
//#define NDEBUG  // disable asserts


/**
 * Custom assert.
 */
#ifndef NDEBUG
#define MCS_ASSERT(cond, msg)                                           \
  do { if (!(cond)) {                                                   \
      std::cerr << "MCS: " << __FILE__ << ": line " << __LINE__ << ": " \
                << "Assertion `"#cond"' failed: " << msg << std::endl   \
                << "Aborted" << std::endl;                              \
      std::exit(EXIT_FAILURE);                                          \
    } } while (0)
#else
#define MCS_ASSERT(cond, mgs) do { } while (0)
#endif


/**
 * Fortran interface.
 */
#ifndef MCS_F77_NAME
#define MCS_F77_NAME(name) name ## _
#endif
#ifndef MCS_F77_CALL
#define MCS_F77_CALL(name) name ## _
#endif



/**
 * @brief MCS namespace.
 */
namespace mcs
{
}


#endif
