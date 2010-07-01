/**
 * @file criterion.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_CRITERION_HH
#define MCS_SUBSET_CRITERION_HH


#include <cmath>

#include "../core/util/math.hh"

#include "linear_model.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value>
    class log_lik
    {

    private:

      const Value m_;

      const Value log_m_;

      const Value log_2pi_;


    public:

      log_lik(const linear_model<Value>& lm) :
	m_(lm.get_observation_count()),
	log_m_(std::log(m_)),
	log_2pi_(std::log(2 * core::util::math::constants<Value>::pi))
      {
      }


      Value
      value(const Value rss) const
      {
	return -0.5 * m_ * (log_2pi_ - log_m_ + std::log(rss) + 1);
      }
	

    };


    class bic
    {

    public:

      template<typename Value>
      class instance
      {

      private:

	const Value m_;

	const Value log_m_;

	const Value log_2pi_;


      public:

	instance(const linear_model<Value>& lm) :
	  m_(lm.get_observation_count()),
	  log_m_(std::log(m_)),
	  log_2pi_(std::log(2 * core::util::math::constants<Value>::pi))
	{
	}


	Value
	value(const Value n,
	      const Value rss) const
	{
          const Value p = n + 1;

	  return log_m_ * p + m_ * (log_2pi_ - log_m_ + std::log(rss));
	}

      };


    public:

      template<typename Value>
      instance<Value>
      for_model(const linear_model<Value>& lm) const
      {
	return instance<Value>(lm);
      }

    };


    class aic
    {

    public:

      template<typename Value>
      class instance
      {

      private:

	const Value m_;

	const Value log_m_;

	const Value log_2pi_;


      public:

	instance(const linear_model<Value>& lm) :
	  m_(lm.get_observation_count()),
	  log_m_(std::log(m_)),
	  log_2pi_(std::log(2 * core::util::math::constants<Value>::pi))
	{
	}


	Value
	value(const Value n,
	      const Value rss) const
	{
          const Value p = n + 1;

	  return 2 * p + m_ * (log_2pi_ + std::log(rss) - log_m_ + 1);
	}

      };


    public:

      template<typename Value>
      instance<Value>
      for_model(const linear_model<Value>& lm) const
      {
	return instance<Value>(lm);
      }

    };


    class aicc
    {

    public:

      template<typename Value>
      class instance
      {

      private:

	const aic::instance<Value> aic_;

	const Value m_;


      public:

	instance(const linear_model<Value>& lm) :
	  aic_(lm),
	  m_(lm.get_observation_count())
	{
	}


	Value
	value(const Value n,
	      const Value rss) const
	{
	  const Value p = n + 1;

	  return aic_.value(n, rss) + (2 * p * (p + 1)) / (m_ - p - 1);
	}

      };


    public:

      template<typename Value>
      instance<Value>
      for_model(const linear_model<Value>& lm) const
      {
	return instance<Value>(lm);
      }

    };


    class aicc2
    {

    public:

      template<typename Value>
      class instance
      {

      private:

	const Value m_;


      public:

	instance(const linear_model<Value>& lm) :
	  m_(lm.get_observation_count())
	{
	}


	Value
	value(const Value n,
	      const Value rss) const
	{
	  const Value p = n + 1;

	  return std::log(rss / m_) + (m_ + p) / (m_ - p - 2);
	}

      };


    public:

      template<typename Value>
      instance<Value>
      for_model(const linear_model<Value>& lm) const
      {
	return instance<Value>(lm);
      }

    };


    struct aicu
    {

    public:

      template<typename Value>
      class instance
      {

      private:

	const Value m_;


      public:

	instance(const linear_model<Value>& lm) :
	  m_(lm.get_observation_count())
	{
	}


	Value
	value(const Value n,
	      const Value rss) const
	{
	  const Value p = n + 1;
	
	  return std::log(rss / (m_ - p)) + (m_ + p) / (m_ - p - 2);
	}

      };


    public:

      template<typename Value>
      instance<Value>
      for_model(const linear_model<Value>& lm) const
      {
	return instance<Value>(lm);
      }

    };


    class cp
    {

    public:

      template<typename Value>
      class instance
      {

      private:

	const Value m_;

	const Value rms_;


      public:

	instance(const linear_model<Value>& lm) :
	  m_(lm.get_observation_count()),
	  rms_(lm.rms())
	{
	}


	Value
	value(const Value n,
	      const Value rss) const
	{
	  const Value p = n + 1;
	
	  return rss / rms_ - m_ + 2 * p;
	}

      };


    public:

      template<typename Value>
      instance<Value>
      for_model(const linear_model<Value>& lm) const
      {
	return instance<Value>(lm);
      }

    };


    template<typename Penalty>
    class maic
    {

    public:

      template<typename Value>
      class instance
      {

      private:

	const Value penalty_;

	const log_lik<Value> logL_;


      public:

	instance(Value penalty, const linear_model<Value>& lm) :
	  penalty_(penalty),
	  logL_(lm)
	{
	}


	Value
	value(const Value n,
	      const Value rss) const
	{
          const Value p = n + 1;

	  return penalty_ * n - 2 * logL_.value(rss);
	}

      };


    private:

      Penalty penalty_;


    public:

      maic(Penalty penalty) :
	penalty_(penalty)
      {
      }


      Penalty
      get_penalty() const
      {
	return penalty_;
      }


      template<typename Value>
      instance<Value>
      for_model(const linear_model<Value>& lm) const
      {
	return instance<Value>(Value(penalty_), lm);
      }

    };


  }

}


#endif
