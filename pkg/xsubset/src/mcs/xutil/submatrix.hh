#ifndef _MCS_XUTIL_SUBMATRIX_HH_
#define _MCS_XUTIL_SUBMATRIX_HH_


namespace MCS
{

  namespace xutil
  {


    template<typename DataType>
    class matrix;


    template<typename DataType>
    class __submatrix_pool;



    template<typename DataType>
    class submatrix : public matrix<DataType>
    {
      friend class matrix<DataType>;
      friend class __submatrix_pool<DataType>;

    private:
      typedef DataType   data_type;
      typedef data_type* data_ptr;

      typedef matrix<DataType>   matrix_type;
      typedef matrix_type&       matrix_ref;
      typedef const matrix_type& const_matrix_ref;

      typedef submatrix<DataType>   submatrix_type;
      typedef submatrix_type&       submatrix_ref;
      typedef const submatrix_type& const_submatrix_ref;

    private:
      submatrix();
      submatrix(data_ptr base, int ldim, int nrow, int ncol);
      submatrix(const_submatrix_ref s);
      ~submatrix();

      void* operator new(size_t);

    public:
      submatrix_ref operator =(data_type x);
      submatrix_ref operator =(const_matrix_ref m);
      submatrix_ref operator =(const_submatrix_ref m);

      void swap(matrix_ref m);
    };


  }

}


#include "submatrix.cc"
#endif
