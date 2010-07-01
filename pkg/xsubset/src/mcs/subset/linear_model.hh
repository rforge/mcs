/**
 * @file linear_model.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_LINEAR_MODEL_HH
#define MCS_SUBSET_LINEAR_MODEL_HH


#include "../core/numeric/matrix.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value>
    class linear_model
    {

    public:

      typedef Value value_type;

      typedef core::numeric::matrix<Value> matrix_type;

      typedef typename matrix_type::size_type size_type;


    private:

      typedef typename matrix_type::range_type range_type;


    private:

      const size_type observation_count_;

      const size_type regressor_count_;

      matrix_type rz_;


    public:

      linear_model(const matrix_type& ay);

      size_type
      get_observation_count() const;

      size_type
      get_regressor_count() const;

      matrix_type
      get_rz() const;

      value_type
      rss() const;

      value_type
      rms() const;

    };


  }

}


#include "linear_model.cc"
#endif
