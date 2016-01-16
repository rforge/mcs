#ifndef R_LM_SUBSETS_CC
#define R_LM_SUBSETS_CC


#include "mcs/subset/dca.hh"
#include "mcs/subset/bba.hh"
#include "mcs/subset/pbba.hh"
#include "mcs/subset/hbba.hh"
#include "mcs/subset/xbba.hh"
#include "mcs/subset/criteria.hh"

#include <cstring>


extern "C"
void
R_lmSubsets(const char* const* const algo, const int* const nobs,
            const int* const nvar, const int* const size,
            const int* const mark, const int* const nbest,
            const int* const pmin, const int* const v,
            const double* const xy, int* sDf, double* sRss, int* sWhich,
            const double* const tau, int* const info, int* const nodes)
{
  using namespace mcs::subset;


  // workspace
  int    wIndex [*size * *nbest              ];
  double wRss   [*size * *nbest              ];
  int    wSubset[*size * *nbest * (*size + 1)];


  *info = 0;

  if (std::strcmp(*algo, "dca") == 0)
    {
      dca(*nobs, *size, *mark, *nbest, v, xy,
          *nobs, wIndex, wRss, wSubset, *nodes);
    }
  else if (std::strcmp(*algo, "bba") == 0)
    {
      bba(*nobs, *size, *mark, *nbest, v, xy,
          *nobs, wIndex, wRss, wSubset, *nodes);
    }
  else if (std::strcmp(*algo, "pbba") == 0)
    {
      pbba(*nobs, *size, *mark, *nbest, *pmin, v, xy,
           *nobs, wIndex, wRss, wSubset, *nodes);
    }
  else if (std::strcmp(*algo, "hbba") == 0)
    {
      hbba(*nobs, *size, *mark, *nbest, *pmin, v, xy,
           *nobs, wIndex, wRss, wSubset, tau, *nodes);
    }
  else if (std::strncmp(*algo, "xbba", 4) == 0)
    {
      xbba(*algo, *nobs, *size, *mark, *nbest, *pmin, v, xy,
           *nobs, wIndex, wRss, wSubset, *info, *nodes);
    }
  else
    {
      *info = -1;
    }


  if (*info >= 0)
    {
      for (int i = *mark; i < *size; ++i)
        {
          for (int j = 0; j < *nbest; ++j)
            {
              int offset = i * *nbest + j;
              int index  = wIndex[offset];

	      const int n = wSubset[index * (*size + 1)];

	      if (n > 0) {
		sDf [offset] = n + 1;
		sRss[offset] = wRss[index];

		offset *= *nvar;
		index  *= (*size + 1);
		for (int k = 0; k < i + 1; ++k)
		  {
		    int v = wSubset[index + k + 1];

		    sWhich[offset + v] = 1;
		  }
	      }
	    }
        }
    }
}



extern "C"
void
R_lmSelect(const char* const* const algo, const int* const nobs,
           const int* const nvar, const int* const size,
           const int* const mark, const int* const nbest,
           const int* const pmin, const int* const v,
           const double* const xy, int* sDf, double* sRss, double* sVal,
	   int* sWhich, const double* penalty, const double* const tau,
           int* const info, int* const nodes)
{
  using namespace mcs::subset;


  Criteria::Aic<double> aic(*penalty, *nobs);


  // workspace
  int    wIndex [*nbest              ];
  double wRss   [*nbest              ];
  double wVal   [*nbest              ];
  int    wSubset[*nbest * (*size + 1)];


  *info = 0;

  if (std::strcmp(*algo, "dca") == 0)
    {
      dca(*nobs, *size, *mark, *nbest, v, xy,
          *nobs, wIndex, wRss, wVal, wSubset, aic, *nodes);
    }
  else if (std::strcmp(*algo, "bba") == 0)
    {
      bba(*nobs, *size, *mark, *nbest, v, xy,
          *nobs, wIndex, wRss, wVal, wSubset, aic, *nodes);
    }
  else if (std::strcmp(*algo, "pbba") == 0)
    {
      pbba(*nobs, *size, *mark, *nbest, *pmin, v, xy,
           *nobs, wIndex, wRss, wVal, wSubset, aic, *nodes);
    }
  else if (std::strcmp(*algo, "hbba") == 0)
    {
      hbba(*nobs, *size, *mark, *nbest, *pmin, v, xy,
           *nobs, wIndex, wRss, wVal, wSubset, aic, *tau, *nodes);
    }
  else if (std::strncmp(*algo, "xbba", 4) == 0)
    {
      xbba(*algo, *nobs, *size, *mark, *nbest, *pmin, v, xy,
           *nobs, wIndex, wRss, wVal, wSubset, aic, *info, *nodes);
    }
  else
    {
      *info = -1;
    }


  if (*info >= 0)
    {
      for (int i = 0; i < *nbest; ++i)
        {
          int offset = i;
          int index  = wIndex[offset];

	  const int n = wSubset[index * (*size + 1)];

	  if (n > 0) {
	    sDf [offset] = n + 1;
	    sRss[offset] = wRss[index];
	    sVal[offset] = wVal[index];

	    offset *= *nvar;
	    index  *= (*size + 1);
	    for (int j = 0; j < wSubset[index]; ++j)
	      {
		int v = wSubset[index + j + 1];

		sWhich[offset + v] = 1;
	      }
	  }
        }
    }
}



#endif
