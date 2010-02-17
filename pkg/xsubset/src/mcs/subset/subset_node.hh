/**
 * @file subset_node.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_SUBSET_NODE_HH
#define MCS_SUBSET_SUBSET_NODE_HH


#include <vector>

#include "../core/numeric/givens.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Subset,
	     typename Matrix>
    class subset_node
    {

    public:

      typedef Subset subset_type;

      typedef typename subset_type::variable_type variable_type;

      typedef typename subset_type::value_type value_type;

      typedef typename subset_type::size_type size_type;

      typedef Matrix matrix_type;

      typedef typename Subset::const_iterator const_subset_iterator;


    private:

      typedef typename Matrix::range_type Range;

      typedef core::numeric::givens<value_type> Givens;


    private:

      subset_type subset_;

      matrix_type rz_;

      size_type mark_;


      std::vector<size_type> work_sizes_;
      std::vector<value_type> work_values_;
      std::vector<Givens> work_givens_;
      matrix_type work_matrix_;


    public:

      subset_node(size_type max_subset_size);

      subset_node(const subset_type& s,
		  const matrix_type& rz,
		  size_type mark);

      size_type
      get_subset_size() const;

      size_type
      get_mark() const;

      value_type
      get_rss() const;

      const_subset_iterator
      subset_begin() const;

      const_subset_iterator
      subset_end() const;

      value_type
      get_z(size_type i) const;
      
      void
      drop_regressor(size_type j,
		     subset_node& res) const;

      void
      preorder_regressors();

      void
      swap(subset_node& other);


    private:

      void
      init_work();

    };


  }

}


#include "subset_node.cc"
#endif
