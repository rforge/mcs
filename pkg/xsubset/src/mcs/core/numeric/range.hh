/**
 * @file range.hh
 *
 * @author Marc Hofmann
 */

#ifndef _MCS_CORE_NUMERIC_RANGE_HH_
#define _MCS_CORE_NUMERIC_RANGE_HH_


#include <limits>

#include "../../mcs.hh"


#define MAX_LEN std::numeric_limits<size_type>::max()


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Size = size_t>
      class range
      {

      public:

        typedef Size size_type;


      private:

        size_type pos_;

        size_type len_;


      public:

	range(size_type begin) :
	  pos_(begin),
	  len_(MAX_LEN)
	{
	  MCS_ASSERT(begin >= 0);
	}

        range(size_type begin,
              size_type end) :
          pos_(begin),
          len_(end - begin)
        {
          MCS_ASSERT(begin >= 0);
          MCS_ASSERT(end >= begin);
        }

	bool
	open() const
	{
	  return len_ == MAX_LEN;
	}

        size_type
        pos() const
        {
          return pos_;
        }

        size_type
        len() const
        {
          return len_;
        }

      };


    }

  }

}


#endif
