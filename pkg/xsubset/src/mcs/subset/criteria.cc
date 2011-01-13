/**
 * @file critera.cc
 */
#ifndef MCS_SUBSET_CRITERA_HH
#define MCS_SUBSET_CRITERA_HH


#include <cmath>

#include "../core/util/math.hh"

#include "lm.hh"
#include "criteria.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    log_lik<Value, Size>::log_lik(const Size nobs) :
      nobs_(nobs),
      log_nobs_(std::log(nobs_)),
      log_2pi_(std::log(2 * core::util::math::constants<Value>::pi))
    {
    }


    template<typename Value,
             typename Size>
    Value
    log_lik<Value, Size>::value(const Value rss) const
    {
      return -0.5 * nobs_ * (log_2pi_ - log_nobs_ + std::log(rss) + 1);
    }



    template<typename Value,
             typename Size>
    aic<Value, Size>::instance::instance(const Value K,
                                         const Size nobs) :
      K_(K),
      L_(nobs)
    {
    }


    template<typename Value,
             typename Size>
    Value
    aic<Value, Size>::instance::value(const Size nvar,
                                      const Value rss) const
    {
      // nvar + 1  to account for sd as estimated parameter
      const Size npar = nvar + 1;

      return -2 * L_.value(rss) + K_ * npar;
    }


    template<typename Value,
             typename Size>
    aic<Value, Size>::aic() :
      K_(2)
    {
    }


    template<typename Value,
             typename Size>
    aic<Value, Size>::aic(const Value K) :
      K_(K)
    {
    }


    template<typename Value,
             typename Size>
    typename aic<Value, Size>::instance
    aic<Value, Size>::for_model(const lm<Value, Size>& x) const
    {
      return instance(K_, x.nobs());
    }



    template<typename Value,
             typename Size>
    cp<Value, Size>::instance::instance(const Value M,
                                        const Size nobs,
                                        const Value rms) :
      M_(M),
      nobs_(nobs),
      rms_(rms)
    {
    }


    template<typename Value,
             typename Size>
    Value
    cp<Value, Size>::instance::value(const Size nvar,
                                     const Value rss) const
    {
      return rss / rms_ - (nobs_ - M_ * nvar);
    }


    template<typename Value,
             typename Size>
    cp<Value, Size>::cp() :
      M_(2)
    {
    }


    template<typename Value,
             typename Size>
    cp<Value, Size>::cp(const Value M) :
      M_(M)
    {
    }


    template<typename Value,
             typename Size>
    typename cp<Value, Size>::instance
    cp<Value, Size>::for_model(const lm<Value, Size>& x) const
    {
      return instance(M_, x.nobs(), x.rms());
    }


  }

}


#endif
