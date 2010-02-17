/**
 * @file crit_selector.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_CRIT_SELECTOR_HH
#define MCS_SUBSET_CRIT_SELECTOR_HH


#include <vector>

#include "../core/numeric/matrix.hh"

#include "subset.hh"
#include "subset_node.hh"
#include "selector_base.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value,
	     typename Crit>
    class crit_selector;


    template<typename Value,
	     typename Crit>
    struct selector_traits<crit_selector<Value, Crit> >
    {

      typedef Value value_type;

    };


    template<typename Value,
	     typename Crit>
    class crit_selector :
      public selector_base<crit_selector<Value, Crit> >
    {

      friend class selector_base<crit_selector<Value, Crit> >;


    private:

      typedef crit_selector<Value, Crit> Self;

      typedef selector_base<Self> Base;


    public:

      typedef typename Base::value_type value_type;

      typedef typename Base::matrix_type matrix_type;

      typedef typename Base::size_type size_type;

      typedef typename Base::subset_type subset_type;

      typedef typename Base::linear_model_type linear_model_type;

      typedef typename Crit::template instance<value_type> criterion_type;


    protected:

      using Base::current_node_;

      value_type tolerance_;

      subset_type best_subset_;

      const criterion_type crit_;


    public:

      crit_selector(const linear_model_type& lm,
		    size_type mark);

      crit_selector(const linear_model_type& lm,
		    size_type mark,
		    const Crit& crit);

      void
      set_tolerance(value_type tolerance);

      subset_type
      get_best_subset() const;


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


#include "crit_selector.cc"
#endif
