/**
 * @file lm.hh
 */
#ifndef MCS_SUBSET_LM_HH
#define MCS_SUBSET_LM_HH


#include "../core/numeric/matrix.hh"


namespace mcs
{

  namespace subset
  {


    template<typename Value,
             typename Size>
    class lm
    {

    private:

      const Size nobs_;

      const Size nvar_;

      mcs::core::numeric::matrix<Value, Size> rz_;

      Value rss_;


    public:

      lm(Size nobs,
         Size nvar,
         const Value* ay);

      lm(const mcs::core::numeric::matrix<Value, Size>& ay);

      Size
      nobs() const;

      Size
      nvar() const;

      mcs::core::numeric::matrix<Value, Size>
      rz() const;

      Value
      rss() const;

      Value
      rms() const;

    };


  }

}


#include "lm.cc"
#endif




