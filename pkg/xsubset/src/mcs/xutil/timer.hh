#ifndef _MCS_XROBUST_TIMER_HH_
#define _MCS_XROBUST_TIMER_HH_


#include <ctime>
#include <limits>


namespace MCS
{

  namespace xrobust
  {


    // BOOST
    class timer
    {
    private:
      std::clock_t _start_time;

    public:
      timer() {
        _start_time = std::clock();
      }

      void restart()
      {
        _start_time = std::clock();
      }

      double elapsed() const
      {
        return  double(std::clock() - _start_time) / CLOCKS_PER_SEC;
      }

      double elapsed_max() const
      {
        return (double((std::numeric_limits<std::clock_t>::max)())
                - double(_start_time)) / double(CLOCKS_PER_SEC);
      }

      double elapsed_min() const
      {
        return double(1)/double(CLOCKS_PER_SEC);
      }

    };


  }

}


#endif
