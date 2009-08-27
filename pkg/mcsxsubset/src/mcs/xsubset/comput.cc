#ifndef _MCS_XSUBSET_COMPUT_CC_
#define _MCS_XSUBSET_COMPUT_CC_


#include <cmath>
#include <limits>
#include <algorithm>

#include "../xutil/assert.hh"
#include "../xutil/blas.hh"
#include "../xutil/lapack.hh"
#include "../xutil/sort.hh"

#include "xsubset.hh"
#include "comput.hh"


namespace MCS
{

  namespace xsubset
  {


    using namespace xutil;


    double DBL_MAX = std::numeric_limits<double>::max();


    double
    sqr(double dx)
    {
      return (dx * dx);
    }


    double
    sign(double dx)
    {
      if (dx > 0)
        {
          return 1.0L;
        }
      else if (dx < 0)
        {
          return -1.0L;
        }
      else
        {
          return 0.0L;
        }
    }


    double
    copysign(double dx, double dy)
    {
      return (dx * sign(dy));
    }


    t_rot
    rot_gen(double dx, double dy)
    {
      t_rot g;

      if (dy == 0.0L)
        {
          g.c = copysign(1.0L, dx);
          g.s = 0.0L;
          g.r = std::abs(dx);
        }
      else if (dx == 0)
        {
          g.c = 0.0L;
          g.s = copysign(1.0L, dy);
          g.r = std::abs(dy);
        }
      else if (std::abs(dy) > std::abs(dx))
        {
          double t = dx / dy;
          double u = copysign(std::sqrt(1 + t * t), dy);
          g.s = 1 / u;
          g.c = g.s * t;
          g.r = dy * u;
        }
      else
        {
          double t = dy / dx;
          double u = copysign(std::sqrt(1 + t * t), dx);
          g.c = 1 / u;
          g.s = g.c * t;
          g.r = dx * u;
        }

      return g;
    }


    void
    chol_init(t_chol& c, t_size size)
    {
      c.size = size;
      c.mark = 0;
      c.rz = t_dmatrix(size, size, 0.0L);
    }


    void
    chol_init(t_chol& c, const t_dmatrix& ay)
    {
      c.size = ay.ncol();
      c.mark = 0;
      c.rz = ay;
      lapack::geqrf(c.rz);
    }


    void
    chol_copy(const t_chol& c1, t_chol& c2)
    {
      c2.size = c1.size;
      c2.mark = c1.mark;
      lapack::lacpy("Upper",
                    range(c1.mark, c1.size),
                    range(c1.mark, c1.size),
                    c1.rz, c2.rz);
    }


    double
    chol_rss(const t_chol& c)
    {
      return sqr(c.rz.diag()(c.size - 1));
    }


    void
    chol_rss(const t_chol& c, t_dbl_iter rss)
    {
      const t_size size = c.size;
      const t_size mark = c.mark;
      const t_size n = size - mark - 1;

      rss[n - 1] = sqr(c.rz(size - 1,
                            size - 1));
      for (t_size j = 2; j <= n; ++j)
        {
          rss[n - j] =
            rss[n - j + 1]      +
            sqr(c.rz(size - j,
                     size - 1));
        }
    }


    void
    chol_swap(t_chol& c1, t_chol& c2)
    {
      std::swap(c1.size, c2.size);
      std::swap(c1.mark, c2.mark);
      c1.rz.swap(c2.rz);
    }


    void
    chol_shift(t_chol& c)
    {
      ++(c.mark);
    }


    void
    chol_shift(t_chol& c, t_size n)
    {
      c.mark += n;
    }


    void
    chol_shift(const t_chol& c1, t_chol& c2)
    {
      chol_copy(c1, c2);
      chol_shift(c2);
    }


    void
    chol_drop(t_chol& c)
    {
      --(c.size);

      t_size i = c.mark;

      t_rot g = rot_gen(c.rz(i    , i + 1),
                        c.rz(i + 1, i + 1));
      c.rz(i    , i) = g.r;
      c.rz(i + 1, i) = 0.0L;

      for (t_size j = i + 1; j < c.size; ++j)
        {
          c.rz(i, j) =
            g.c * c.rz(i    , j + 1) +
            g.s * c.rz(i + 1, j + 1);
          c.rz(i + 1, j) =
            -g.s * c.rz(i    , j + 1) +
            g.c  * c.rz(i + 1, j + 1);
        }

      for (i = c.mark + 1; i < c.size; ++i)
        {
          g = rot_gen(c.rz(i    , i    ),
                      c.rz(i + 1, i + 1));
          c.rz(i    , i) = g.r;
          c.rz(i + 1, i) = 0.0L;

          for (t_size j = i + 1; j < c.size; ++j)
            {
              double t =
                g.c * c.rz(i    , j    ) +
                g.s * c.rz(i + 1, j + 1);
              c.rz(i + 1, j) =
                -g.s * c.rz(i    , j    ) +
                g.c  * c.rz(i + 1, j + 1);
              c.rz(i, j) = t;
            }
        }
    }


    void
    chol_drop(const t_chol& c1, t_chol& c2)
    {
      c2.size = c1.size - 1;
      c2.mark = c1.mark;

      t_size i = c2.mark;

      t_rot g = rot_gen(c1.rz(i    , i + 1),
                        c1.rz(i + 1, i + 1));
      c2.rz(i    , i) = g.r;
      c2.rz(i + 1, i) = 0.0L;

      for (t_size j = i + 1; j < c2.size; ++j)
        {
          c2.rz(i, j) =
            g.c * c1.rz(i    , j + 1) +
            g.s * c1.rz(i + 1, j + 1);
          c2.rz(i + 1, j) =
            -g.s * c1.rz(i    , j + 1) +
            g.c  * c1.rz(i + 1, j + 1);
        }

      for (i = c2.mark + 1; i < c2.size; ++i)
        {
          g = rot_gen(c2.rz(i    ,     i),
                      c1.rz(i + 1, i + 1));
          c2.rz(i    , i) = g.r;
          c2.rz(i + 1, i) = 0.0L;

          for (t_size j = i + 1; j < c2.size; ++j)
            {
              double t =
                g.c * c2.rz(i    ,     j) +
                g.s * c1.rz(i + 1, j + 1);
              c2.rz(i + 1, j) =
                -g.s * c2.rz(i    , j    ) +
                g.c  * c1.rz(i + 1, j + 1);
              c2.rz(i, j) = t;
            }
        }
    }


    void
    chol_bounds(const t_chol& c, t_dbl_iter rss)
    {
      t_rot g[c.size];

      for (t_size i = c.mark; i < (c.size - 1); ++i)
        {
          for (t_size j = i + 1; j < c.size; ++j)
            {
              double t = c.rz(i, j);
              for (t_size k = i + 1; k < j; ++k)
                {
                  t =
                    -g[k].s * t          +
                    g[k].c  * c.rz(k, j);
                }
              g[j] = rot_gen(t, c.rz(j, j));
            }

          rss[i - c.mark] = sqr(g[c.size - 1].r);
        }
    }


    void
    chol_permute(t_chol& c, t_size_iter index, t_dmatrix& work)
    {
      const t_size size = c.size;
      const t_size mark = c.mark;

      work.swap(c.rz);
      for (t_size i = mark; i < (size - 1); ++i)
        {
          t_size j = mark + index[i - mark];
          blas::copy(range(mark, j + 1),
                     work.col(j),
                     c.rz.col(i));
          c.rz(range(j + 1, size), i) = 0.0L;
        }
      blas::copy(range(mark, size),
                 work.col(size - 1),
                 c.rz.col(size - 1));
      lapack::geqrf(range(mark, size),
                    range(mark, size),
                    c.rz);
    }


  }

}


#endif
