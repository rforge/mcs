/**
 * @file crit_selector.cc
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_CRIT_SELECTOR_CC
#define MCS_SUBSET_CRIT_SELECTOR_CC


#include <cmath>
#include <algorithm>

#include "crit_selector.hh"


#define CRIT_SELECTOR crit_selector<Value, Crit>


namespace mcs
{

  namespace subset
  {
    

    template<typename Value,
	     typename Crit>
    CRIT_SELECTOR::crit_selector(const linear_model_type& lm,
				 const size_type mark) :
      Base::selector_base(lm, mark),
      tolerance_(0),
      best_subset_(lm.get_regressor_count()),
      crit_(lm)
    {
    }
    

    template<typename Value,
	     typename Crit>
    CRIT_SELECTOR::crit_selector(const linear_model_type& lm,
				 const size_type mark,
				 const Crit& crit) :
      Base::selector_base(lm, mark),
      tolerance_(0),
      best_subset_(lm.get_regressor_count()),
      crit_(crit.for_model(lm))
    {
    }


    template<typename Value,
	     typename Crit>
    void
    CRIT_SELECTOR::set_tolerance(const value_type tolerance)
    {
      tolerance_ = tolerance;
    }


    template<typename Value,
	     typename Crit>
    typename CRIT_SELECTOR::subset_type
    CRIT_SELECTOR::get_best_subset() const
    {
      return best_subset_;
    }


    template<typename Value,
	     typename Crit>
    void
    CRIT_SELECTOR::init()
    {
      Base::init_impl();
    }


    template<typename Value,
	     typename Crit>
    void
    CRIT_SELECTOR::next_node()
    {
      Base::next_node_impl();
    }


    template<typename Value,
	     typename Crit>
    void
    CRIT_SELECTOR::preorder_regressors()
    {
      Base::preorder_regressors_impl();
    }


    template<typename Value,
	     typename Crit>
    void
    CRIT_SELECTOR::report_subsets()
    {
      namespace math = core::util::math;

      const size_type n = current_node_.get_subset_size();
      const size_type k = current_node_.get_mark();

      const auto first = current_node_.subset_begin();
      auto last = current_node_.subset_end();

      value_type rss = current_node_.get_rss();
      for (size_type j = n - 1; j >= k; --j, --last)
	{
	  const value_type val = crit_.value(j + 1, rss);
	  if (val < best_subset_.get_value())
	    {
	      best_subset_.assign(first, last);
	      best_subset_.set_value(val);
	    }
	  rss += math::sqr(current_node_.get_z(j));
	}
    }


    template<typename Value,
	     typename Crit>
    void
    CRIT_SELECTOR::spawn_nodes()
    {
      Base::spawn_nodes_impl();
    }


    template<typename Value,
	     typename Crit>
    bool
    CRIT_SELECTOR::cut_child(const size_type k) const
    {
      const value_type rss = current_node_.get_rss();

      const value_type val = crit_.value(k + 1, rss);
      const value_type bound = (1 + tolerance_) * val;

      return bound >= best_subset_.get_value();
    }


  }

}


#undef CRIT_SELECTOR


#endif
