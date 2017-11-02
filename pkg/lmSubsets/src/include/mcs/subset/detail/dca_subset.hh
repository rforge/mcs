#ifndef MCS_SUBSET_DETAIL_DCA_SUBSET_HH
#define MCS_SUBSET_DETAIL_DCA_SUBSET_HH



#include <algorithm>  // std::copy_n, std::rotate
#include <iterator>  // std::back_inserter, std::make_reverse_iterator
#include <utility>  // std::swap
#include <vector>



#include "mcs/util/algo.hh"  // arrange



namespace mcs    {
namespace subset {
namespace detail {



class dca_subset
{

public:

    static void
    drop_column(
        const std::vector<int>& subset,
        const int mark,
        std::vector<int>& result
    ) noexcept
    {
        const int n = subset.size();
        const int k = mark;

        const auto first = subset.begin();
        const auto skip = first + k;
        const auto last = subset.end();

        result.assign(first, skip);
        result.insert(result.end(), skip + 1, last);
    }



    static void
    permute_complete(
        const std::vector<int>& subset,
        const int mark,
        const std::vector<int>& pos,
        std::vector<int>& result
    ) noexcept
    {
        using namespace mcs;

        const int n = subset.size();
        const int k = mark;
        const int p = n - k;

        result.assign(subset.begin(), subset.begin() + k);
        util::arrange_n(subset.begin() + k, p, pos.begin(),
                        std::back_inserter(result));
    }



    static void
    permute_partial_1(
        std::vector<int>& subset,
        const int mark,
        const int pos
    ) noexcept
    {
        using std::swap;

        const int k = mark;
        const int j = pos;

        swap(subset[k], subset[k + j]);
    }



    static void
    permute_partial_2(
        std::vector<int>& subset,
        const int mark,
        const int pos
    ) noexcept
    {
        const int k = mark;
        const int j = pos;

        std::rotate(std::make_reverse_iterator(subset.begin() + k + j + 1),
                    std::make_reverse_iterator(subset.begin() + k + j),
                    std::make_reverse_iterator(subset.begin() + k));
    }



public:

    dca_subset() = delete;

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
