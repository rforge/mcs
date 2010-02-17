/**
 * @file subset_node.cc
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_SUBSET_NODE_CC
#define MCS_SUBSET_SUBSET_NODE_CC


#include <cmath>
#include <algorithm>

#include "../core/numeric/blas.hh"
#include "../core/numeric/lapack.hh"
#include "../core/numeric/givens.hh"

#include "../core/util/algo.hh"
#include "../core/util/math.hh"

#include "subset_node.hh"


#define SUBSET_NODE subset_node<Subset, Matrix>


namespace mcs
{

  namespace subset
  {


    template<typename Subset,
	     typename Matrix>
    SUBSET_NODE::subset_node(const size_type max_subset_size) :
      subset_(max_subset_size),
      rz_(max_subset_size + 1, max_subset_size + 1, value_type(0)),
      mark_()
    {
      init_work();
    }


    template<typename Subset,
	     typename Matrix>
    SUBSET_NODE::subset_node(const subset_type& s,
			     const matrix_type& rz,
			     const size_type mark) :
      subset_(s),
      rz_(rz),
      mark_(mark)
    {
    }


    template<typename Subset,
	     typename Matrix>
    typename SUBSET_NODE::size_type
    SUBSET_NODE::get_subset_size() const
    {
      return subset_.get_size();
    }


    template<typename Subset,
	     typename Matrix>
    typename SUBSET_NODE::size_type
    SUBSET_NODE::get_mark() const
    {
      return mark_;
    }


    template<typename Subset,
	     typename Matrix>
    typename SUBSET_NODE::const_subset_iterator
    SUBSET_NODE::subset_begin() const
    {
      return subset_.begin();
    }


    template<typename Subset,
	     typename Matrix>
    typename SUBSET_NODE::const_subset_iterator
    SUBSET_NODE::subset_end() const
    {
      return subset_.end();
    }


    template<typename Subset,
	     typename Matrix>
    typename SUBSET_NODE::value_type
    SUBSET_NODE::get_z(const size_type i) const
    {
      const size_type n = get_subset_size();

      return rz_(i, n);
    }


    template<typename Subset,
	     typename Matrix>
    typename SUBSET_NODE::value_type
    SUBSET_NODE::get_rss() const
    {
      return subset_.get_value();
    }


    template<typename Subset,
	     typename Matrix>
    void
    SUBSET_NODE::drop_regressor(const size_type j,
				subset_node& res) const
    {
      namespace math = core::util::math;

      const size_type n = get_subset_size();

      const auto first = subset_.begin();
      const auto drop = first + j;
      const auto last = subset_.end();
      res.subset_.assign(first, drop);
      res.subset_.append(drop + 1, last);

      Givens::zero(    rz_(j    , Range(j + 1, n + 1)),
		       rz_(j + 1, Range(j + 1, n + 1)),
		   res.rz_(j    , Range(j    , n    )),
		   res.rz_(j + 1, Range(j    , n    )));

      for (size_type i = j + 1; i < n; ++i)
	{
	  Givens::zero(res.rz_(i    , Range(i    , n    )),
		           rz_(i + 1, Range(i + 1, n + 1)),
		       res.rz_(i    , Range(i    , n    )),
		       res.rz_(i + 1, Range(i    , n    )));
	}

      res.mark_ = j;

      res.subset_.set_value(math::sqr(res.rz_(n - 1, n - 1)));
    }


    template<typename Subset,
	     typename Matrix>
    void
    SUBSET_NODE::preorder_regressors()
    {
      namespace algo = core::util::algo;
      namespace blas = core::numeric::blas;
      namespace lapack = core::numeric::lapack;

      const size_type n = get_subset_size();
      const size_type k = get_mark();
      const size_type d = n - k;

      for (size_type j = k; j < n; ++j)
      	{
      	  for (size_type i = j + 1; i < n + 1 ; ++i)
      	    {
      	      value_type t = rz_(j, i);
      	      for (size_type g = j + 1; g < i; ++g)
      		{
      		  t = -work_givens_[g].s() * t + work_givens_[g].c() * rz_(g, i);
      		}
      	      work_givens_[i].gen(t, rz_(i, i));
      	    }

      	  work_values_[j - k] = std::abs(work_givens_[n].r());
      	}

      algo::iota_n(work_sizes_.begin(), d, 0);
      algo::index_sort_n(work_sizes_.begin(), d,
			 work_values_.begin(), std::greater<value_type>());

      algo::permute_n(subset_.begin() + k, d, work_sizes_.begin());

      work_matrix_.swap(rz_);
      for (size_type i = k; i < n; ++i)
      	{
      	  size_type j = k + work_sizes_[i - k];
      	  blas::copy(Range(k, j + 1), work_matrix_.col(j), rz_.col(i));
      	  rz_(Range(j + 1, n + 1), i).fill(0);
      	}
      blas::copy(Range(k, n + 1), work_matrix_.col(n), rz_.col(n));
      lapack::geqrf(Range(k, n + 1), rz_);
    }


    template<typename Subset,
	     typename Matrix>
    void
    SUBSET_NODE::swap(subset_node& other)
    {
      std::swap(subset_, other.subset_);
      std::swap(rz_, other.rz_);
      std::swap(mark_, other.mark_);
    }


    template<typename Subset,
	     typename Matrix>
    void
    SUBSET_NODE::init_work()
    {
      const size_type nn = subset_.get_max_size() + 1;

      work_sizes_.resize(nn);
      work_values_.resize(nn);
      work_givens_.resize(nn);
      work_matrix_ = matrix_type(nn, nn);
    }


  }

}


#undef SUBSET_NODE


#endif
