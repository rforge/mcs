#ifndef _MCS_XUTIL_IO_CC_
#define _MCS_XUTIL_IO_CC__


#include <iostream>
#include <cmath>

#include "io.hh"
#include "vector.hh"
#include "matrix.hh"


namespace MCS
{

  namespace xutil
  {



    template<typename T>
    std::istream&
    operator >>(std::istream& is, vector<T>& v)
    {
      for (int i = 0; i < v.len(); ++i)
	{
          is >> v(i);
	}
      return is;
    }


    template<typename T>
    std::ostream&
    operator <<(std::ostream& os, const vector<T>& v)
    {
      for (int i = 0; i < v.len(); ++i)
        {
          os << " " << v(i);
        }
      os << std::endl;
      return os;
    }


    template<typename T>
    std::ostream&
    operator <<(std::ostream& os, const matrix<T>& m)
    {
      int i, j;

      for (i = 0; i < m.nrow(); ++i)
        {
          for (j = 0; j < m.ncol(); ++j)
            {
              os << " " << m(i, j);
            }
          os << std::endl;
        }
      return os;
    }


    template<typename T>
    std::istream&
    operator >>(std::istream& is, matrix<T>& m)
    {
      int i, j;

      for (i = 0; i < m.nrow(); ++i)
        {
          for (j = 0; j < m.ncol(); ++j)
            {
              is >> m(i, j);
            }
         }
      return is;
    }


  }

}


#endif
