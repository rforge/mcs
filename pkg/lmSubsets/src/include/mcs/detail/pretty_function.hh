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
