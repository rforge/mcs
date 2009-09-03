#ifndef _MCS_XUTIL_MATRIX_HH_
#define _MCS_XUTIL_MATRIX_HH_


namespace MCS
{

  namespace xutil
  {


    class range;


    template<typename DataType>
    class vector;


    template<typename DataType>
    class subvector;


    template<typename DataType>
    class submatrix;



    template<typename DataType>
    class matrix
    {
      friend class submatrix<DataType>;

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

      typedef matrix<data_type>  matrix_type;
      typedef matrix_type&       matrix_ref;
      typedef const matrix_type& const_matrix_ref;

      typedef submatrix<data_type>  submatrix_type;
      typedef submatrix_type&       submatrix_ref;
      typedef const submatrix_type& const_submatrix_ref;

    private:
      data_ptr _stor;

      data_ptr _base;
      int _ldim;
      int _nrow;
      int _ncol;

    public:
      matrix();
      matrix(int nrow, int ncol);
      matrix(int nrow, int ncol, data_type x);
      matrix(int nrow, int ncol, data_ptr x);
      matrix(const_matrix_ref m);

    private:
      matrix(data_ptr stor, int nrow, int ncol);
      matrix(data_ptr base, int ldim, int nrow, int ncol);

    public:
      ~matrix();

      matrix_ref operator =(data_type x);
      matrix_ref operator =(matrix_type m);

      data_ref operator ()(int i, int j);
      const_data_ref operator ()(int i, int j) const;

      subvector_ref operator ()(int i, range rj);
      const_subvector_ref operator ()(int i, range rj) const;

      subvector_ref operator ()(range ri, int j);
      const_subvector_ref operator ()(range ri, int j) const;

      submatrix_ref operator ()(range ri, range rj);
      const_submatrix_ref operator ()(range ri, range rj) const;

      int ldim() const;
      int nrow() const;
      int ncol() const;

      subvector_ref row(int i);
      const_subvector_ref row(int i) const;

      submatrix_ref row(range ri);
      const_submatrix_ref row(range ri) const;

      subvector_ref col(int j);
      const_subvector_ref col(int j) const;

      submatrix_ref col(range rj);
      const_submatrix_ref col(range rj) const;

      subvector_ref diag();
      const_subvector_ref diag() const;

      void swap(matrix_ref m);
      void swap(submatrix_ref m);

    private:
      void copy(data_type x);
      void copy(const_matrix_ref m);

      void move(matrix_ref m);
    };


  }

}


#include "matrix.cc"
#endif
