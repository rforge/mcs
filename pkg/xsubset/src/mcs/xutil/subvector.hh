#ifndef _MCS_XUTIL_SUBVECTOR_HH_
#define _MCS_XUTIL_SUBVECTOR_HH_


namespace MCS
{

  namespace xutil
  {


    template<typename DataType>
    class vector;


    template<typename DataType>
    class matrix;


    template<typename DataType>
    class __subvector_pool;



    template<typename DataType>
    class subvector : public vector<DataType>
    {
      friend class vector<DataType>;
      friend class matrix<DataType>;

      friend class __subvector_pool<DataType>;

    private:
      typedef DataType   data_type;
      typedef data_type* data_ptr;

      typedef vector<DataType>   vector_type;
      typedef vector_type&       vector_ref;
      typedef const vector_type& const_vector_ref;

      typedef subvector<DataType>   subvector_type;
      typedef subvector_type&       subvector_ref;
      typedef const subvector_type& const_subvector_ref;

    private:
      subvector();
      subvector(data_ptr base, int inc, int len);
      subvector(const_subvector_ref s);
      ~subvector();

      void* operator new(size_t);

    public:
      subvector_ref operator =(data_type x);
      subvector_ref operator =(const_vector_ref v);
      subvector_ref operator =(const_subvector_ref v);

      void swap(vector_ref v);
    };


  }

}


#include "subvector.cc"
#endif
