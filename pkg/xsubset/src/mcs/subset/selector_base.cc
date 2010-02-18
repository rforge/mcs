/**
 * @file selector_base.cc
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_SELECTOR_BASE_CC
#define MCS_SUBSET_SELECTOR_BASE_CC


#include "selector_base.hh"


#define SELECTOR_BASE selector_base<Selector>


namespace mcs
{

  namespace subset
  {


    template<typename Selector>
    SELECTOR_BASE::selector_base(const linear_model_type& lm,
				 const size_type mark) :
      preordering_radius_(0),
      lm_(lm),
      mark_(mark),
      node_container_(lm.get_regressor_count(),
		      Node(lm.get_regressor_count())),
      current_node_(lm.get_regressor_count())
    {
    }


    template<typename Selector>
    void
    SELECTOR_BASE::set_preordering_radius(const size_type preordering_radius)
    {
      preordering_radius_ = preordering_radius;
    }


    template<typename Selector>
    unsigned int
    SELECTOR_BASE::get_node_count() const
    {
      return node_count_;
    }


    template<typename Selector>
    void
    SELECTOR_BASE::select()
    {
      init();

      while (node_stack_top_ != node_stack_bottom_)
	{

	  ++node_count_;

	  next_node();

	  preorder_regressors();
	  report_subsets();

	  spawn_nodes();

	}
    }


    template<typename Selector>
    void
    SELECTOR_BASE::init()
    {
      static_cast<Selector*>(this)->init();
    }


    template<typename Selector>
    void
    SELECTOR_BASE::next_node()
    {
      static_cast<Selector*>(this)->next_node();
    }


    template<typename Selector>
    void
    SELECTOR_BASE::preorder_regressors()
    {
      const size_type n = current_node_.get_subset_size();
      const size_type k = current_node_.get_mark();
      const size_type d = n - k;
      
      if (d > preordering_minimum_d_)
	{
	  static_cast<Selector*>(this)->preorder_regressors();
	}
    }


    template<typename Selector>
    void
    SELECTOR_BASE::report_subsets()
    {
      static_cast<Selector*>(this)->report_subsets();
    }


    template<typename Selector>
    void
    SELECTOR_BASE::spawn_nodes()
    {
      static_cast<Selector*>(this)->spawn_nodes();
    }


    template<typename Selector>
    bool
    SELECTOR_BASE::cut_child(const size_type k) const
    {
      return static_cast<const Selector*>(this)->cut_child(k);
    }
      

    template<typename Selector>
    void
    SELECTOR_BASE::init_impl()
    {
      const size_type n = lm_.get_regressor_count();
      const size_type k = mark_;
      const size_type d = n - k;

      if (preordering_radius_ > d)
	{
	  preordering_minimum_d_ = 0;
	}
      else
	{
	  preordering_minimum_d_ = d - preordering_radius_;
	}

      node_stack_bottom_ = node_container_.begin();
      node_stack_top_ = node_stack_bottom_;

      subset_type s(n, lm_.rss());

      *node_stack_top_ = Node(s, lm_.get_rz(), k);
      ++node_stack_top_;

      node_count_ = 0;
    }
      

    template<typename Selector>
    void
    SELECTOR_BASE::next_node_impl()
    {
      current_node_.swap(*(--node_stack_top_));
    }


    template<typename Selector>
    void
    SELECTOR_BASE::preorder_regressors_impl()
    {
      const size_type n = current_node_.get_subset_size();
      const size_type k = current_node_.get_mark();
      const size_type d = n - k;

      if (d > preordering_minimum_d_)
	{
	  current_node_.preorder_regressors();
	}
    }


    template<typename Selector>
    void
    SELECTOR_BASE::report_subsets_impl()
    {
      // subclass repsonsibility
    }


    template<typename Selector>
    void
    SELECTOR_BASE::spawn_nodes_impl()
    {
      const size_type n = current_node_.get_subset_size();
      const size_type k = current_node_.get_mark();

      for (size_type j = k; j < (n - 1); ++j)
	{
	  if (!cut_child(j))
	    {
	      current_node_.drop_regressor(j, *(node_stack_top_++));
	    }
	}
    }


    template<typename Selector>
    bool
    SELECTOR_BASE::cut_child_impl(const size_type k) const
    {
      // subclass responsibility
    }


  }

}


#undef SELECTOR_BASE


#endif
