/**
 * @file criteria.hh
 */
#ifndef MCS_SUBSET_CRITERIA_HH
#define MCS_SUBSET_CRITERIA_HH


#include "lm.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    class log_lik
    {

    private:

      const Value nobs_;

      const Value log_nobs_;

      const Value log_2pi_;


    public:

      log_lik(Size nobs);

      Value
      value(Value rss) const;

    };



    template<typename Value,
             typename Size>
    class aic
    {

    public:

      class instance
      {

      private:

        const Value K_;

	const log_lik<Value, Size> L_;


      public:

	instance(Value K,
                 Size nobs);

	Value
	value(const Size nvar,
	      const Value rss) const;

      };


    private:

      const Value K_;


    public:

      aic();

      aic(Value K);

      instance
      for_model(const lm<Value, Size>& x) const;

    };



    template<typename Value,
             typename Size>
    class cp
    {

    public:

      class instance
      {

      private:

        const Value M_;

	const Value nobs_;

	const Value rms_;


      public:

	instance(Value M_,
                 Size nobs,
                 Value rms);

	Value
	value(Size nvar,
	      Value rss) const;

      };


    private:

      const Value M_;


    public:

      cp();

      cp(Value M);

      instance
      for_model(const lm<Value, Size>& x) const;

    };


  }

}


#include "criteria.cc"
#endif
