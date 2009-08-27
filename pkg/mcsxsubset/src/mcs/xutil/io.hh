#ifndef _MCS_XUTIL_IO_HH_
#define _MCS_XUTIL_IO_HH_


#include <iostream>

#include "vector.hh"
#include "matrix.hh"


namespace MCS
{

  namespace xutil
  {


    template<typename T>
    std::istream&
    operator >>(std::istream& is, vector<T>& v);

    template<typename T>
    std::ostream&
    operator <<(std::ostream& os, const vector<T>& v);



    template<typename T>
    std::ostream&
    operator <<(std::ostream& os, const matrix<T>& m);

    template<typename T>
    std::istream&
    operator >>(std::istream& is, matrix<T>& m);


  }

}


#include "io.cc"
#endif
