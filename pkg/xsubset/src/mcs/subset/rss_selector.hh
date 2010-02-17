/**
 * @file rss_selector.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_RSS_SELECTOR_HH
#define MCS_SUBSET_RSS_SELECTOR_HH


#include <vector>

#include "selector_base.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value>
    class rss_selector;


    template<typename Value>
    struct selector_traits<rss_selector<Value> >
    {

      typedef Value value_type;

    };


    template<typename Value>
    class rss_selector :
      public selector_base<rss_selector<Value> >
    {

      friend class selector_base<rss_selector<Value> >;


    private:

      typedef rss_selector<Value> Self;

      typedef selector_base<Self> Base;


    public:

      typedef typename Base::value_type value_type;

      typedef typename Base::matrix_type matrix_type;

      typedef typename Base::size_type size_type;

      typedef typename Base::subset_type subset_type;

      typedef typename Base::linear_model_type linear_model_type;


    private:

      typedef std::vector<value_type> ValueContainer;

      typedef std::vector<subset_type> SubsetContainer;


    protected:

      using Base::current_node_;

      ValueContainer tolerances_;

      SubsetContainer best_subsets_;

      size_type cutting_mark_;


    public:

      rss_selector(const linear_model_type& lm,
		   size_type mark);

      template<typename ValueIterator>
      ValueIterator
      load_tolerances(ValueIterator src);

      template<typename SubsetIterator>
      SubsetIterator
      store_best_subsets(SubsetIterator dst) const;

      subset_type
      get_best_subset(size_type subset_size) const;


    protected:

      void
      init();

      void
      next_node();

      void
      preorder_regressors();

      void
      report_subsets();

      void
      spawn_nodes();

      bool
      cut_child(size_type k) const;

    };


  }

}


#include "rss_selector.cc"
#endif
