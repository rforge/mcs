// Copyright 2018  Marc Hofmann
//
// This file is part of 'lmSubsets'.
//
// 'lmSubsets' is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// 'lmSubsets' is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with 'lmSubsets'.  If not, see <http://www.gnu.org/licenses/>.



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
