#ifndef _MCS_XUTIL_BLAS_HH_
#define _MCS_XUTIL_BLAS_HH_


#include "assert.hh"

#include "vector.hh"
#include "matrix.hh"


extern "C"
{


  void
  dcopy_(const int* n,
	 const double* x, int* incx,
	 const double* y, int* incy);


  void
  dgemv_(const char* trans,
         const int* m, const int* n,
         const double* alpha,
         const double* a, const int* lda,
         const double* x, const int* incx,
         const double* beta,
         double* y, const int* incy);


  void
  drot_(const int* n,
	double* x, const int* incx,
	double* y, const int* incy,
	const double* c, const double* s);


  void
  drotg_(double* x, double* y,
	 double* c, double *s);

}


namespace MCS
{

  namespace xutil
  {

    namespace blas
    {


      void
      copy(int n,
           const double& x, int incx,
           double& y, int incy)
      {
        dcopy_(&n, &x, &incx, &y, &incy);
      }


      void
      copy(int n,
           const vector<double>& x,
           vector<double>& y)
      {
        MCS_ASSERT((n <= x.len()) && (n <= y.len()));

        copy(n,
             x(0), x.inc(),
             y(0), y.inc());
      }


      void
      copy(const vector<double>& x,
           vector<double>& y)
      {
        MCS_ASSERT(x.len() == y.len());

        copy(x.len(), x, y);
      }


      void
      copy(const range& r,
           const vector<double>& x,
           vector<double>& y)
      {
        copy(x(r), y(r));
      }


      void
      gemv(const char* trans, int m, int n,
           double alpha,
           const double& a, int lda,
           const double& x, int incx,
           double beta,
           double& y, int incy)
      {
        dgemv_(trans, &m, &n, &alpha, &a, &lda,
               &x, &incx, &beta, &y, &incy);
      }


      void
      gemv(const char* trans,
           double alpha,
           const matrix<double>& a,
           const vector<double>& x,
           double beta,
           vector<double>& y)
      {
        MCS_ASSERT(a.ncol() == x.len());
        MCS_ASSERT(a.nrow() == y.len());

        gemv(trans, a.nrow(), a.ncol(), alpha,
             a(0, 0), a.ldim(), x(0), x.inc(),
             beta, y(0), y.inc());
      }


      void
      rot(int n,
          double& x, int incx,
          double& y, int incy,
          double c, double s)
      {
        drot_(&n, &x, &incx, &y, &incy, &c, &s);
      }


      void
      rot(int n,
          vector<double>& x,
          vector<double>& y,
          double c, double s)
      {
        MCS_ASSERT(x.len() >= n);
        MCS_ASSERT(y.len() >= n);

        rot(n, x(0), x.inc(), y(0), y.inc(), c, s);
      }


      void
      rot(vector<double>& x,
          vector<double>& y,
          double c, double s)
      {
        MCS_ASSERT(x.len() == y.len());

        rot(x.len(), x(0), x.inc(), y(0), y.inc(), c, s);
      }


      void
      rotg(double& x, double& y, double& c, double& s)
      {
        drotg_(&x, &y, &c, &s);
      }


    }

  }

}


#endif
