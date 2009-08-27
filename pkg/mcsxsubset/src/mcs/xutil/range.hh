#ifndef _MCS_XUTIL_RANGE_HH_
#define _MCS_XUTIL_RANGE_HH_


#include "assert.hh"


namespace MCS
{

  namespace xutil
  {


    template<typename T>
    class vector;


    template<typename T>
    class matrix;


    class range
    {
      template<typename T>
      friend class vector;

      template<typename T>
      friend class matrix;

    private:
      int _start;
      int _end;

    public:
      range(int start, int end = -1) : _start(start), _end(end)
      {
        MCS_ASSERT(start >= 0);
        MCS_ASSERT(end == -1 || end >= start);
      }
    };


  }

}


#endif
