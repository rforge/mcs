#ifndef R_MCS_SUBSET_CC
#define R_MCS_SUBSET_CC


#include "mcs/subset/dca.hh"
#include "mcs/subset/bba.hh"
#include "mcs/subset/pbba.hh"
#include "mcs/subset/hbba.hh"
#include "mcs/subset/xbba.hh"
#include "mcs/subset/criteria.hh"

#include <cstring>



extern "C"
void
R_mcsSubset(const char* const* const algo, const int* const nobs,
            const int* const nvar, const int* const size,
            const int* const mark, const int* const nbest,
            const int* const pmin, const int* const v,
            const double* const xy, double* rss, double* aic, int* which,
            const double* penalty, const double* const tau,
            int* const nodes, int* const info)
{
  using namespace mcs::subset;


  if (*penalty == 0)
    {
      int    sIndex[(*size - *mark) * *nbest        ];
      double sRss  [(*size - *mark) * *nbest        ];
      int    s     [(*size - *mark) * *nbest * *size];

      *info = 0;

      if (std::strcmp(*algo, "dca") == 0)
        {
          *nodes = dca(*nobs, *size, *mark, *nbest, v, xy, *nobs, sIndex, sRss, s);
        }
      else if (std::strcmp(*algo, "bba") == 0)
        {
          *nodes = bba(*nobs, *size, *mark, *nbest, v, xy, *nobs, sIndex, sRss, s);
        }
      else if (std::strcmp(*algo, "pbba") == 0)
        {
          *nodes = pbba(*nobs, *size, *mark, *nbest, *pmin, v, xy, *nobs, sIndex, sRss, s);
        }
      else if (std::strcmp(*algo, "hbba") == 0)
        {
          *nodes = hbba(*nobs, *size, *mark, *nbest, *pmin, v, xy, *nobs, sIndex, sRss, s, tau + *mark);
        }
      else if (std::strncmp(*algo, "xbba", 4) == 0)
        {
          *nodes = xbba(*algo, info, *nobs, *size, *mark, *nbest, *pmin, v, xy, *nobs, sIndex, sRss, s);
        }
      else
        {
          *info = 1;
        }

      if (*info == 0)
        {
          for (int i = *mark; i < *size; ++i)
            {
              for (int j = 0; j < *nbest; ++j)
                {
                  int offset = (i - *mark) * *nbest + j;
                  int index  = sIndex[offset];

                  offset += *mark * *nbest;

                  rss[offset] = sRss[index];

                  offset *= *nvar;
                  index  *= *size;
                  for (int k = 0; k < i + 1; ++k)
                    {
                      int v = s[index + k];

                      if (v >= 0) which[offset + v] = 1;
                    }
                }
            }
        }
    }
  else
    {
      Criteria::Aic<double> c(*penalty, *nobs);

      int    sIndex[*nbest        ];
      double sRss  [*nbest        ];
      double sCrit [*nbest        ];
      int    sSize [*nbest        ];
      int    s     [*nbest * *size];

      *info = 0;

      if (std::strcmp(*algo, "dca") == 0)
        {
          *nodes = dca(*nobs, *size, *mark, *nbest, v, xy, *nobs, sIndex, sRss, sCrit, sSize, s, c);
        }
      else if (std::strcmp(*algo, "bba") == 0)
        {
          *nodes = bba(*nobs, *size, *mark, *nbest, v, xy, *nobs, sIndex, sRss, sCrit, sSize, s, c);
        }
      else if (std::strcmp(*algo, "pbba") == 0)
        {
          *nodes = pbba(*nobs, *size, *mark, *nbest, *pmin, v, xy, *nobs, sIndex, sRss, sCrit, sSize, s, c);
        }
      else if (std::strcmp(*algo, "hbba") == 0)
        {
          *nodes = hbba(*nobs, *size, *mark, *nbest, *pmin, v, xy, *nobs, sIndex, sRss, sCrit, sSize, s, c, *tau);
        }
      else if (std::strncmp(*algo, "xbba", 4) == 0)
        {
          *nodes = xbba(*algo, info, *nobs, *size, *mark, *nbest, *pmin, v, xy, *nobs, sIndex, sRss, sCrit, sSize, s, c);
        }
      else
        {
          *info = 1;
        }

      if (*info == 0)
        {
          for (int i = 0; i < *nbest; ++i)
            {
              int offset = i;
              int index  = sIndex[offset];
              int n      = sSize[index];

              rss[offset] = sRss[index];
              aic[offset] = sCrit[index];

              offset *= *nvar;
              index  *= *size;
              for (int j = 0; j < n; ++j)
                {
                  int v = s[index + j];

                  if (v >= 0) which[offset + v] = 1;
                }
            }
        }
    }
}



#endif
