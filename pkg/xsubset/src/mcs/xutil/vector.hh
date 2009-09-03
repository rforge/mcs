#ifndef _MCS_XUTIL_VECTOR_HH_
#define _MCS_XUTIL_VECTOR_HH_


#include <algorithm>


namespace MCS
{

  namespace xutil
  {


    class range;


    template<typename DataType>
    class matrix;


    template<typename DataType>
    class subvector;



    template<typename DataType>
    class vector
    {
      friend class matrix<DataType>;
      friend class subvector<DataType>;

    private:
      typedef DataType         data_type;
      typedef data_type*       data_ptr;
      typedef data_type&       data_ref;
      typedef const data_type& const_data_ref;

      typedef vector<data_type>  vector_type;
      typedef vector_type&       vector_ref;
      typedef const vector_type& const_vector_ref;

      typedef subvector<data_type>  subvector_type;
      typedef subvector_type&       subvector_ref;
      typedef const subvector_type& const_subvector_ref;

    private:
      data_ptr _stor;

      data_ptr _base;
      int _inc;
      int _len;

    public:
      vector();
      vector(int len);
      vector(int len, data_type x);
      vector(const_vector_ref v);

    private:
      vector(data_ptr stor, int len);
      vector(data_ptr base, int inc, int len);

    public:
      ~vector();

      vector_ref operator =(data_type x);
      vector_ref operator =(vector_type v);

      data_ref operator ()(int i);
      const_data_ref operator ()(int i) const;

      subvector_ref operator ()(range r);
      const_subvector_ref operator ()(range r) const;

      int inc() const;
      int len() const;

      void swap(vector_ref v);
      void swap(subvector_ref v);

    private:
      void copy(data_type x);
      void copy(const_vector_ref v);

      void move(vector_ref v);
    };


  }

}


#include "vector.cc"
#endif
