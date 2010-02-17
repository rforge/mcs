/**
 * @file matrix.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_NUMERIC_MATRIX_HH
#define MCS_CORE_NUMERIC_MATRIX_HH


#include <memory>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Size>
      class range;


      template<typename Value,
	       typename Alloc>
      class vector;


      template<typename Value,
	       typename Alloc>
      class subvector;


      template<typename Value,
	       typename Alloc>
      class submatrix;


      template<typename Value,
	       typename Alloc = std::allocator<Value> >
      class matrix
      {

	friend class submatrix<Value, Alloc>;


      public:

	typedef Value value_type;

        typedef Alloc allocator_type;

        typedef typename allocator_type::size_type size_type;

        typedef typename allocator_type::difference_type difference_type;

        typedef typename allocator_type::pointer pointer;

        typedef typename allocator_type::const_pointer const_pointer;

        typedef typename allocator_type::reference reference;

        typedef typename allocator_type::const_reference const_reference;

        typedef range<size_type> range_type;

        typedef subvector<value_type, allocator_type> subvector_type;

        typedef submatrix<value_type, allocator_type> submatrix_type;


      private:

	struct matrix_impl : public allocator_type
	{

	  pointer stor;

	  pointer base;

	  size_type ldim;

	  size_type nrow;

	  size_type ncol;


	  matrix_impl();

	  matrix_impl(const allocator_type& alloc);

	};


	matrix_impl impl_;


      public:

	matrix();

	matrix(size_type nrow,
	       size_type ncol);

	template<typename V>
	matrix(size_type nrow,
	       size_type ncol,
	       V val);

	template<typename V>
	matrix(size_type nrow,
	       size_type ncol,
	       const V* data);

	matrix(const matrix& mat);

	matrix(matrix&& mat);

	matrix(submatrix_type&& mat);


      private:

	matrix(pointer stor,
	       size_type nrow,
	       size_type ncol);

	matrix(pointer base,
	       size_type ldim,
	       size_type nrow,
	       size_type ncol);


      public:

	~matrix();

	matrix&
	operator =(matrix mat);

	reference
	operator ()(size_type irow,
		    size_type icol);

	const_reference
	operator ()(size_type irow,
		    size_type icol) const;

	subvector_type
	operator ()(size_type irow,
		    range_type rcol);

	const subvector_type
	operator ()(size_type irow,
		    range_type rcol) const;

	subvector_type
	operator ()(range_type rrow,
		    size_type icol);

	const subvector_type
	operator ()(range_type rrow,
		    size_type icol) const;

	submatrix_type
	operator ()(range_type rrow,
		    range_type rcol);

	const submatrix_type
	operator ()(range_type rrow,
		    range_type rcol) const;

	size_type
	ldim() const;

	size_type
	nrow() const;

	size_type
	ncol() const;

	subvector_type
	row(size_type irow);

	const subvector_type
	row(size_type irow) const;

	subvector_type
	col(size_type icol);

	const subvector_type
	col(size_type icol) const;

	submatrix_type
	row(range_type rrow);

	const submatrix_type
	row(range_type rrow) const;

	submatrix_type
	col(range_type rcol);

	const submatrix_type
	col(range_type rcol) const;

	subvector_type
	diag();

	const subvector_type
	diag() const;

	template<typename V>
	void
	fill(V val);

	template<typename V>
	void
	fill(const V* val);

	void
	copy(const matrix& mat);

	void
	swap(matrix& mat);


      private:

	reference
	at(size_type irow,
	   size_type icol) const;

      };


      template<typename Value,
	       typename Alloc>
      class submatrix : public matrix<Value, Alloc>
      {

	friend class matrix<Value, Alloc>;


      public:

	typedef Value value_type;

        typedef Alloc allocator_type;

        typedef typename allocator_type::size_type size_type;

        typedef typename allocator_type::difference_type difference_type;

        typedef typename allocator_type::pointer pointer;

        typedef typename allocator_type::const_pointer const_pointer;

        typedef typename allocator_type::reference reference;

        typedef typename allocator_type::const_reference const_reference;

        typedef range<size_type> range_type;

        typedef subvector<value_type, allocator_type> subvector_type;

        typedef matrix<value_type, allocator_type> matrix_type;


      private:

	submatrix(pointer base,
		  size_type ldim,
		  size_type nrow,
		  size_type ncol);


      public:

	~submatrix();

	submatrix&
	operator =(const matrix_type& mat);

	submatrix&
	operator =(matrix_type&& mat);

      };


    }

  }

}


#define MATRIX mcs::core::numeric::matrix<Value, Alloc>


namespace std
{


  template<typename Value,
	   typename Alloc>
  void
  swap(MATRIX& x, MATRIX& y);


}


#undef MATRIX


#include "matrix.cc"
#endif
