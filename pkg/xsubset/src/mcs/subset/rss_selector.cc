/**
 * @file rss_selector.cc
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_RSS_SELECTOR_CC
#define MCS_SUBSET_RSS_SELECTOR_CC


#include <iostream>
#include <algorithm>

#include "rss_selector.hh"


#define RSS_SELECTOR rss_selector<Value>


namespace mcs
{

  namespace subset
  {
    

    template<typename Value>
    RSS_SELECTOR::rss_selector(const linear_model_type& lm,
			       const size_type mark) :
      Base::selector_base(lm, mark),
      tolerances_(lm.get_regressor_count(), 0),
      best_subsets_(lm.get_regressor_count(),
		    subset_type(lm.get_regressor_count()))
    {
    }


    template<typename Value>
    template<typename ValueIterator>
    ValueIterator
    RSS_SELECTOR::load_tolerances(ValueIterator src)
    {
      std::copy_n(src, tolerances_.size(), tolerances_.begin());
      return src + tolerances_.size();
    }


    template<typename Value>
    template<typename SubsetIterator>
    SubsetIterator
    RSS_SELECTOR::store_best_subsets(SubsetIterator dst) const
    {
      return std::copy(best_subsets_.begin(), best_subsets_.end(), dst);
    }


    template<typename Value>
    typename RSS_SELECTOR::subset_type
    RSS_SELECTOR::get_best_subset(size_type subset_size) const
    {
      return best_subsets_[subset_size - 1];
    }


    template<typename Value>
    void
    RSS_SELECTOR::init()
    {
      Base::init_impl();
    }


    template<typename Value>
    void
    RSS_SELECTOR::next_node()
    {
      Base::next_node_impl();

      const size_type n = current_node_.get_subset_size();
      const size_type k = current_node_.get_mark();

      const value_type rss = current_node_.get_rss();

      for (cutting_mark_ = n - 2;
	   cutting_mark_ >= k;
	   --cutting_mark_)
	{
	  const value_type bound = (1 + tolerances_[cutting_mark_]) * rss;
	  if (bound < best_subsets_[cutting_mark_].get_value())
	    {
	      break;
	    }
	}
    }


    template<typename Value>
    void
    RSS_SELECTOR::preorder_regressors()
    {
      Base::preorder_regressors_impl();
    }


    template<typename Value>
    void
    RSS_SELECTOR::report_subsets()
    {
      namespace math = core::util::math;

      const size_type n = current_node_.get_subset_size();
      const size_type k = current_node_.get_mark();
      const size_type d = n - k;

      value_type rss = current_node_.get_rss();
      const auto first = current_node_.subset_begin();
      auto last = current_node_.subset_end();
      for (size_type  j = n - 1; j >= k; --j, --last)
	{
	  if (rss < best_subsets_[j].get_value())
	    {
	      best_subsets_[j].assign(first, last);
	      best_subsets_[j].set_value(rss);
	    }
	  rss += math::sqr(current_node_.get_z(j));
	}
    }


    template<typename Value>
    void
    RSS_SELECTOR::spawn_nodes()
    {
      Base::spawn_nodes_impl();
    }


    template<typename Value>
    bool
    RSS_SELECTOR::cut_child(const size_type k) const
    {
      return k > cutting_mark_;
    }


  }

}


#undef RSS_SELECTOR


#endif
