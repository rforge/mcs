#ifndef MCS_SUBSET_DETAIL_DCA_QRZ_HH
#define MCS_SUBSET_DETAIL_DCA_QRZ_HH



#include <algorithm>  // std::copy_n, std::fill_n
#include <cstdlib>  // std::abs
#include <utility>  // std::swap
#include <vector>



#include "mcs/core/givens.hh"
#include "mcs/core/lapack.hh"
#include "mcs/core/matrix.hh"



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar>
class dca_qrz
{

    using matrix = mcs::core::matrix<Scalar>;

    using matrix_span = mcs::core::matrix<Scalar&>;

    using matrix_cspan = mcs::core::matrix<const Scalar&>;

    using givens = mcs::core::givens<Scalar>;

    using lapack = mcs::core::lapack<Scalar>;



private:

    mutable std::vector<Scalar> aux_work_;

    mutable std::vector<Scalar> aux_tau_;

    mutable std::vector<givens> aux_givens_;



public:

    dca_qrz(const int root_size) noexcept :
        aux_work_(root_size + 1),
        aux_tau_(root_size + 1),
        aux_givens_(root_size + 1)
    {
    }



public:

    matrix
    rz(matrix_cspan ay_mat) noexcept
    {
        const int n = ay_mat.ncol() - 1;

        matrix rz_tmp = ay_mat;
        lapack::geqr2(rz_tmp, aux_tau_, aux_work_);

        return rz_tmp({0, n + 1}, {0, n + 1});
    }



    void
    drop_column(
        matrix_cspan rz_mat,
        const int mark,
        matrix_span out_mat
    ) const noexcept
    {
        const int n = rz_mat.ncol() - 1;
        const int k = mark;

        drop_column(n - k, rz_mat.ptr(k, k), rz_mat.ldim(), out_mat.ptr(k, k),
                    out_mat.ldim());
    }



    void
    column_bounds(
        matrix_cspan rz_mat,
        const int mark,
        std::vector<Scalar>& out
    ) const noexcept
    {
        const int n = rz_mat.ncol() - 1;
        const int k = mark;
        const int p = n - k;

        column_bounds(p, rz_mat.ptr(k, k), rz_mat.ldim(), out.data(),
                      aux_givens_.data());
    }



    void
    permute_complete(
        matrix_cspan rz_mat,
        const int mark,
        const std::vector<int>& pos,
        matrix_span out_mat
    ) const noexcept
    {
        const int n = rz_mat.ncol() - 1;
        const int k = mark;
        const int p = n - k;

        permute_complete(p, rz_mat.ptr(k, k), rz_mat.ldim(), pos.data(),
                         out_mat.ptr(k, k), out_mat.ldim(), aux_tau_.data(),
                         aux_work_.data());
    }


    void
    permute_partial_1(
        matrix_span rz_mat,
        const int mark,
        const int pos
    ) const noexcept
    {
        const int n = rz_mat.ncol() - 1;
        const int k = mark;

        permute_partial_1(n - k, rz_mat.ptr(k, k), rz_mat.ldim(), pos - k);
    }



    void
    permute_partial_2(
        matrix_span rz_mat,
        const int mark,
        const int pos
    ) const noexcept
    {
        const int n = rz_mat.ncol() - 1;
        const int k = mark;

        permute_partial_2(n - k, rz_mat.ptr(k, k), rz_mat.ldim(), pos - k);
    }



private:

    static void
    drop_column(
        int n,
        const Scalar* rz,
        const int ldrz,
        Scalar* out,
        const int ldout
    ) noexcept
    {
        givens::zero(n, rz + ldrz, ldrz, rz + ldrz + 1, ldrz, out, ldout,
                     out + 1, ldout);

        while (--n > 0)
        {
            rz += ldrz + 1;
            out += ldout + 1;

            givens::zero(n, out, ldout, rz + ldrz + 1, ldrz, out, ldout,
                         out + 1, ldout);
        }
    }



    static void
    column_bounds(
        const int n,
        const Scalar* const rz,
        int ldrz,
        Scalar* out,
        givens* const aux_givens
    ) noexcept
    {
        const Scalar* colj = rz;

        for (int j = 0; j < n; ++j, colj += ldrz)
        {
            const Scalar* coli = colj + ldrz;

            for (int i = j + 1; i <= n; ++i, coli += ldrz)
            {
                Scalar t = coli[j];

                for (int g = j + 1; g < i; ++g)
                {
                    t = -aux_givens[g].s() * t + aux_givens[g].c() * coli[g];
                }

                aux_givens[i].gen(t, coli[i]);
            }

            *(out++) = std::abs(aux_givens[n].r());
        }
    }



    static void
    permute_complete(
        const int n,
        const Scalar* rz,
        const int ldrz,
        const int* const pos,
        Scalar* const out,
        const int ldout,
        Scalar* const aux_tau,
        Scalar* const aux_work
    ) noexcept
    {
        for (int i = 0; i < n; ++i)
        {
            const int j = pos[i];

            std::copy_n(rz + j * ldrz, j + 1, out + i * ldout);
            std::fill_n(out + i * ldout + j + 1, n - j, 0);
        }

        std::copy_n(rz + n * ldrz, n + 1, out + n * ldout);
        lapack::geqr2(n + 1, n + 1, out, ldout, aux_tau, aux_work);
    }



    static void
    permute_partial_1(
        const int n,
        Scalar* const rz,
        const int ldrz,
        const int pos
    ) noexcept
    {
        const int j = pos;

        Scalar* const col0 = rz;
        Scalar* const colj = rz + j * ldrz;

        for (int i = j - 1; i >= 0; --i)
        {
            const givens g(colj[i], colj[i + 1]);
            g.rot(j - i, colj - ldrz + i, -ldrz, colj - ldrz + i + 1, -ldrz);
            g.rot(n - j, colj + ldrz + i, ldrz, colj + ldrz + i + 1, ldrz);

            colj[i] = g.r();
            colj[i + 1] = 0;
        }

        std::swap(col0[0], colj[0]);
        std::swap(col0[1], colj[1]);

        Scalar* coli = col0;

        for (int i = 1; i < j; ++i)
        {
            coli += ldrz;

            const givens g(coli[i], coli[i + 1]);
            g.rot(n - i, coli + ldrz + i, ldrz, coli + ldrz + i + 1, ldrz);

            coli[i] = g.r();
            coli[i + 1] = 0;
        }
    }



    static void
    permute_partial_2(
        const int n,
        Scalar* const rz,
        const int ldrz,
        const int pos
    ) noexcept
    {
        if (pos < 1)
        {
            return;
        }

        const int j = pos;

        Scalar* const col0 = rz;
        Scalar* const colj = rz + j * ldrz;

        *(colj - ldrz + j) = 0;

        givens g(colj[j - 1], colj[j]);
        g.rot(1, colj - ldrz + j - 1, -ldrz, colj - ldrz + j, -ldrz,
              colj + j - 1, -ldrz, colj + j, -ldrz);
        g.rot(n - j, colj + ldrz + j - 1, ldrz, colj + ldrz + j, ldrz);

        Scalar* rii = colj - ldrz + j - 1;
        for (int i = j - 1; i > 0; --i, rii -= ldrz + 1)
        {
            *rii = 0;

            g.gen(colj[i - 1], g.r());
            g.rot(j - i + 1, colj - ldrz + i - 1, -ldrz, colj + i, -ldrz,
                  colj + i - 1, -ldrz, colj + i, -ldrz);
            g.rot(n - j, colj + ldrz + i - 1, ldrz, colj + ldrz + i, ldrz);
        }

        col0[0] = g.r();
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
