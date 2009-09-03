#ifndef _MCS_XUTIL_LAPACK_HH_
#define _MCS_XUTIL_LAPACK_HH_


#include <algorithm>

#include "assert.hh"

#include "vector.hh"
#include "matrix.hh"


extern "C"
{


  void
  dlacpy_(const char* uplo,
          const int* m, const int* n,
	  const double* a, const int* lda,
	  double* b, const int* ldb);

  void
  dlarf_(const char* side,
         const int* m, const int* n,
         const double* v, const int* incv,
         const double* tau,
         double* c, const int* ldc,
         double* work);

  void
  dlarfg_(const int* n, double* alpha,
          double* x, const int* incx,
          double* tau);


  void
  dgeqrf_(const int* m, const int* n,
          double* a, const int* lda,
          double* tau,
          double* work, const int* lwork,
          int* info);


  void
  dtrtrs_(const char* uplo, const char* trans, const char* diag,
          const int* n, const int* nrhs,
          const double* a, const int* lda,
          double* b, const int* ldb,
          int* info);


}


namespace MCS
{

  namespace xutil
  {

    namespace lapack
    {


      void
      lacpy(const char* uplo, int m, int n,
            const double& a, int lda,
            double& b, int ldb)
      {
        MCS_ASSERT(m >= 0);
        MCS_ASSERT(n >= 0);
        MCS_ASSERT(lda >= std::max(1, m));
        MCS_ASSERT(ldb >= std::max(1, m));

        dlacpy_(uplo, &m, &n,
                &a, &lda,
                &b, &ldb);
      }


      void
      lacpy(const char* uplo, int m, int n,
            const matrix<double>& a,
            matrix<double>& b)
      {
        MCS_ASSERT((a.nrow() >= m) && (a.ncol() >= n));
        MCS_ASSERT((b.nrow() >= m) && (b.ncol() >= n));

        lacpy(uplo, m, n,
              a(0, 0), a.ldim(),
              b(0, 0), b.ldim());
      }


      void
      lacpy(const char* uplo,
            const matrix<double>& a,
            matrix<double>& b)
      {
        MCS_ASSERT(a.nrow() == b.nrow());
        MCS_ASSERT(a.ncol() == b.ncol());

        lacpy(uplo, a.nrow(), a.ncol(),
              a(0, 0), a.ldim(),
              b(0, 0), b.ldim());
      }


      void
      lacpy(const char* uplo,
            const range& r1,
            const range& r2,
            const matrix<double>& a,
            matrix<double>& b)
      {
        lacpy(uplo, a(r1, r2), b(r1, r2));
      }


      void
      larf(const char* side,
           int m, int n,
           const double& v, int incv,
           double tau,
           double& c, int ldc,
           double& work)
      {
        MCS_ASSERT((side[0] == 'L') || (side[0] == 'R'));
        MCS_ASSERT(incv != 0);
        MCS_ASSERT(ldc >= std::max(1, m));

        dlarf_(side,
               &m, &n,
               &v, &incv,
               &tau,
               &c, &ldc,
               &work);
      }


      void
      larf(const char* side,
           const vector<double>& v,
           double tau,
           matrix<double>& c)
      {
        const int m = c.nrow();
        const int n = c.ncol();

        MCS_ASSERT(side[0] == 'L'? v.len() == m : true);
        MCS_ASSERT(side[0] == 'R'? v.len() == n : true);

        double work[std::max(m, n)];

        larf(side,
             m, n,
             v(0), v.inc(),
             tau,
             c(0, 0), c.ldim(),
             work[0]);
      }


      void
      larfg(int n, double& alpha,
            double& x, int incx,
            double& tau)
      {
        dlarfg_(&n, &alpha,
                &x, &incx,
                &tau);
      }


      void
      larfg(double& alpha,
            vector<double>& x,
            double& tau)
      {
        larfg(x.len(), alpha,
              x(0), x.inc(),
              tau);
      }


      void
      geqrf(int m, int n,
            double& a, int lda,
            double& tau,
            double& work, int lwork,
            int& info)
      {
        MCS_ASSERT(m >= 0);
        MCS_ASSERT(n >= 0);
        MCS_ASSERT(lda >= std::max(1, m));
        MCS_ASSERT(lwork >= std::max(1, n));

        dgeqrf_(&m, &n,
                &a, &lda,
                &tau,
                &work, &lwork,
                &info);
      }


      void
      geqrf(matrix<double>& a)
      {
        const int m = a.nrow();
        const int n = a.ncol();
        double tau[std::min(m, n)];
        int lwork = std::max(1, n);
        double work[lwork];
        int info;

        geqrf(m, n,
              a(0, 0), a.ldim(),
              tau[0],
              work[0], lwork,
              info);
      }


      void
      geqrf(const range& r1,
            const range& r2,
            matrix<double>& a)
      {
        geqrf(a(r1, r2));
      }


      void
      trtrs(const char* uplo, const char* trans, const char* diag,
            int n, int nrhs,
            const double& a, int lda,
            double& b, int ldb,
            int& info)
      {
        MCS_ASSERT(n >= 0);
        MCS_ASSERT(nrhs >= 0);
        MCS_ASSERT(lda >= std::max(1, n));
        MCS_ASSERT(ldb >= std::max(1, n));

        dtrtrs_(uplo, trans, diag,
                &n, &nrhs,
                &a, &lda,
                &b, &ldb,
                &info);
      }


      void
      trtrs(const char* uplo, const char* trans, const char* diag,
            const matrix<double>& a,
            matrix<double>& b)
      {
        MCS_ASSERT(a.nrow() == b.nrow());

        int info;
        trtrs(uplo, trans, diag,
              a.nrow(), b.ncol(),
              a(0, 0), a.ldim(),
              b(0, 0), b.ldim(),
              info);
      }


      void
      trtrs(const char* uplo, const char* trans, const char* diag,
            const matrix<double>& a,
            vector<double>& b)
      {
        MCS_ASSERT(a.nrow() == b.len());

        int info;
        trtrs(uplo, trans, diag,
              a.nrow(), 1,
              a(0, 0), a.ldim(),
              b(0), b.len(),
              info);
      }


    }

  }

}


#endif
