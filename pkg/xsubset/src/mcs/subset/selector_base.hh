/**
 * @file selector_base.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_SELECTOR_BASE_HH
#define MCS_SUBSET_SELECTOR_BASE_HH


#include <vector>

#include "../core/numeric/matrix.hh"

#include "subset.hh"
#include "subset_node.hh"
#include "linear_model.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Selector>
    struct selector_traits
    {
    };


    template<typename Selector>
    class selector_base
    {

    public:

      typedef typename selector_traits<Selector>::value_type value_type;

      typedef mcs::core::numeric::matrix<value_type> matrix_type;

      typedef typename matrix_type::size_type size_type;

      typedef subset<size_type, value_type> subset_type;

      typedef linear_model<value_type> linear_model_type;


    private:

      typedef subset_node<subset_type, matrix_type> Node;

      typedef std::vector<Node> NodeContainer;

      typedef typename NodeContainer::iterator NodeIterator;


    protected:

      size_type preordering_radius_;

      size_type preordering_minimum_d_;

      unsigned int node_count_;

      linear_model_type lm_;

      size_type mark_;

      Node current_node_;

      NodeContainer node_container_;

      NodeIterator node_stack_bottom_;

      NodeIterator node_stack_top_;


    protected:

      selector_base(const linear_model_type& lm,
		    const size_type mark);


    public:

      void
      set_preordering_radius(size_type preordering_radius);

      unsigned int
      get_node_count() const;

      void
      select();


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
      cut_child(size_type i) const;

      void
      init_impl();

      void
      next_node_impl();

      void
      preorder_regressors_impl();

      void
      report_subsets_impl();

      void
      spawn_nodes_impl();

      bool
      cut_child_impl(size_type k) const;

    };


  }


}


#include "selector_base.cc"
#endif
