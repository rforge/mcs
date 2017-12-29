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



#ifndef MCS_SUBSET_DETAIL_DCA_RESULT_HH
#define MCS_SUBSET_DETAIL_DCA_RESULT_HH



#include <limits>  // std::numeric_limits
#include <utility>  // std::swap
#include <vector>



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar>
class dca_result
{

    friend void
    swap(dca_result& a, dca_result& b) noexcept
    {
        a.swap(b);
    }



    friend auto
    cbegin(const dca_result& r) noexcept
    {
        return r.cbegin();
    }



    friend auto
    cend(const dca_result& r) noexcept
    {
        return r.cend();
    }



public:

    using value_type = int;



private:

    std::vector<int> subset_;

    Scalar key_;



public:

    dca_result() noexcept :
        key_(std::numeric_limits<Scalar>::quiet_NaN())
    {
    }



    dca_result(
        const std::vector<int>& subset,
        const Scalar key
    ) noexcept :
        subset_(subset),
        key_(key)
    {
        if (subset.size() == 0)
        {
            key_ = std::numeric_limits<Scalar>::quiet_NaN();
        }
    }



public:

    void
    swap(dca_result& other) noexcept
    {
        subset_.swap(other.subset_);
        std::swap(key_, other.key_);
    }



    auto
    cbegin() const noexcept
    {
        return subset_.cbegin();
    }



    auto
    cend() const noexcept
    {
        return subset_.cend();
    }



public:

    operator bool() const noexcept
    {
        return size() > 0;
    }



    int&
    operator [](const int i) noexcept
    {
        return subset_[i];
    }



    const int&
    operator [](const int i) const noexcept
    {
        return subset_[i];
    }



public:

    int
    size() const noexcept
    {
        return subset_.size();
    }



    const std::vector<int>&
    subset() const noexcept
    {
        return subset_;
    }



    Scalar
    key() const noexcept
    {
        return key_;
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
