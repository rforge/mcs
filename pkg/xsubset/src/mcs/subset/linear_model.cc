/**
 * @file linear_model.cc
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_LINEAR_MODEL_CC
#define MCS_SUBSET_LINEAR_MODEL_CC


#include "../core/numeric/lapack.hh"
#include "../core/util/math.hh"

#include "linear_model.hh"


#define LINEAR_MODEL linear_model<Value>


namespace mcs
{

  namespace subset
  {


    template<typename Value>
    LINEAR_MODEL::linear_model(const matrix_type& ay) :
      observation_count_(ay.nrow()),
      regressor_count_(ay.ncol() - 1),
      rz_(ay)
    {
      namespace lapack = core::numeric::lapack;

      lapack::geqrf(rz_);
    }


    template<typename Value>
    typename LINEAR_MODEL::size_type
    LINEAR_MODEL::get_observation_count() const
    {
      return observation_count_;
    }


    template<typename Value>
    typename LINEAR_MODEL::size_type
    LINEAR_MODEL::get_regressor_count() const
    {
      return regressor_count_;
    }


    template<typename Value>
    typename LINEAR_MODEL::matrix_type
    LINEAR_MODEL::get_rz() const
    {
      return rz_;
    }


    template<typename Value>
    typename LINEAR_MODEL::value_type
    LINEAR_MODEL::rss() const
    {
      namespace math = core::util::math;

      const size_type n = regressor_count_;

      return math::sqr(rz_(n, n));
    }


    template<typename Value>
    typename LINEAR_MODEL::value_type
    LINEAR_MODEL::rms() const
    {
      const size_type m = observation_count_;
      const size_type n = regressor_count_;

      return rss() / (m - n);
    }


  }

}


#undef LINEAR_MODEL


#endif
