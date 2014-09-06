#ifndef MCS_SUBSET_CRITERIA_HH
#define MCS_SUBSET_CRITERIA_HH


#include "mcs/subset/detail/log_lik.hh"



namespace mcs    {
namespace subset {


  namespace Criteria {


    template<typename TReal>
    class Aic
    {

    private:

      TReal                 K_;
      detail::LogLik<TReal> L_;


    public:

      Aic(const TReal K, const int nobs) :
        K_{K}   ,
        L_{nobs}
      {
      }


      TReal
      value(const int size, const TReal rss) const
      {
        // size + 1  to account for sd as estimated parameter
        const int npar = size + 1;

        return -2 * L_.value(rss) + K_ * npar;
      }

    };


  }


}
}


#endif
