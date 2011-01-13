/**
 * @file matrix_base.hh
 */
#ifndef MCS_CORE_NUMERIC_MATRIX_BASE_HH
#define MCS_CORE_NUMERIC_MATRIX_BASE_HH


#include "slice.hh"


namespace mcs
{

  namespace core
  {

    namespace numeric
    {



      template<typename Value,
               typename Size>
      class slice_vector;



      template<typename Value,
	       typename Size>
      class matrix;



      template<typename Value,
	       typename Size>
      class slice_matrix;



      template<typename Value,
	       typename Size,
               template<typename V,
                        typename S>
               class Derived>
      class matrix_base
      {

        friend class matrix_base<Value, Size, matrix>;
        friend class matrix_base<Value, Size, slice_matrix>;

        friend class matrix<Value, Size>;
        friend class slice_matrix<Value, Size>;


      protected:

	Value* base_;

        Size ldim_;

	Size nrow_;

	Size ncol_;


      protected:

        matrix_base();

	matrix_base(Value* base,
                    Size ldim,
                    Size nrow,
                    Size ncol);

        matrix_base(const matrix_base<Value, Size, Derived>& x);

	~matrix_base();

      public:

        Value*
        base();

        const Value*
        base() const;

        Size
	ldim() const;

	Size
	nrow() const;

	Size
	ncol() const;

        Value&
        at(Size irow,
           Size icol);

        const Value&
        at(Size irow,
           Size icol) const;

        slice_vector<Value, Size>
        at(Size irow,
           const slice<Size>& scol);

        const slice_vector<Value, Size>
        at(Size irow,
           const slice<Size>& scol) const;

        slice_vector<Value, Size>
        at(const slice<Size>& srow,
           Size icol);

        const slice_vector<Value, Size>
        at(const slice<Size>& srow,
           Size icol) const;

        slice_matrix<Value, Size>
        at(const slice<Size>& srow,
           const slice<Size>& scol);

        const slice_matrix<Value, Size>
        at(const slice<Size>& srow,
           const slice<Size>& scol) const;

	slice_vector<Value, Size>
	row(Size irow);

	const slice_vector<Value, Size>
	row(Size irow) const;

	slice_vector<Value, Size>
	col(Size icol);

	const slice_vector<Value, Size>
	col(Size icol) const;

	slice_matrix<Value, Size>
	rows(const slice<Size>& srow);

	const slice_matrix<Value, Size>
	rows(const slice<Size>& srow) const;

	slice_matrix<Value, Size>
	cols(const slice<Size>& scol);

	const slice_matrix<Value, Size>
	cols(const slice<Size>& rcol) const;

	slice_vector<Value, Size>
	diag();

	const slice_vector<Value, Size>
	diag() const;

        template<template<typename V,
                          typename S>
                 class D>
        void
        copy(const matrix_base<Value, Size, D>& x);

        template<template<typename V,
                          typename S>
                 class D>
        void
        swap(matrix_base<Value, Size, D>&& x);

        void
        fill(Value x);

        matrix_base<Value, Size, Derived>&
        operator =(const matrix_base<Value, Size, Derived>& x);

        template<template<typename V,
                          typename S>
                 class D>
        matrix_base<Value, Size, Derived>&
        operator =(const matrix_base<Value, Size, D>& x);

        template<template<typename V,
                          typename S>
                 class D>
        matrix_base<Value, Size, Derived>&
        operator =(matrix_base<Value, Size, D>&& x);

        matrix_base<Value, Size, Derived>&
        operator =(Value x);

        Value&
	operator ()(Size irow,
		    Size icol);

        const Value&
	operator ()(Size irow,
		    Size icol) const;

	slice_vector<Value, Size>
	operator ()(Size irow,
		    const slice<Size>& scol);

	const slice_vector<Value, Size>
	operator ()(Size irow,
		    const slice<Size>& scol) const;

	slice_vector<Value, Size>
	operator ()(const slice<Size>& srow,
		    Size icol);

	const slice_vector<Value, Size>
	operator ()(const slice<Size>& srow,
		    Size icol) const;

        slice_matrix<Value, Size>
	operator ()(const slice<Size>& srow,
		    const slice<Size>& scol);

	const slice_matrix<Value, Size>
	operator ()(const slice<Size>& srow,
		    const slice<Size>& scol) const;

      };


    }

  }

}


namespace std
{


  template<typename Value,
           typename Size,
           template<typename V,
                    typename S>
           class Derived>
  void
  swap(mcs::core::numeric::matrix_base<Value, Size, Derived>& a,
       mcs::core::numeric::matrix_base<Value, Size, Derived>& b);


}


#include "matrix_base.cc"
#endif
