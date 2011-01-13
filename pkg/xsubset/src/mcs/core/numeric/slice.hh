/**
 * @file slice.hh
 */
#ifndef MCS_CORE_NUMERIC_SLICE_HH
#define MCS_CORE_NUMERIC_SLICE_HH


#include "../../mcs.hh"


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Size>
      class slice
      {

        template<typename V,
                 typename S>
        friend class vector;

        template<typename V,
                 typename S>
        friend class slice_vector;

        template<typename V,
                 typename S,
                 template<typename, typename>
                 class D>
        friend class matrix_base;


      private:

        Size pos_;

        Size len_;


      public:

        slice(Size begin,
              Size end) :
          pos_(begin),
          len_(end - begin)
        {
          MCS_ASSERT(begin >= 0, "invalid argument");
          MCS_ASSERT(end >= begin, "invalid argument");
        }

        Size
        pos() const
        {
          return pos_;
        }

        Size
        len() const
        {
          return len_;
        }

      };


    }

  }

}


#endif
