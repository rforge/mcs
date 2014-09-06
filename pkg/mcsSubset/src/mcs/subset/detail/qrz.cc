#ifndef MCS_SUBSET_DETAIL_QRZ_CC
#define MCS_SUBSET_DETAIL_QRZ_CC


namespace mcs    {
namespace subset {
namespace detail {



  template<typename TReal>
  void
  Qrz<TReal>::drop(int n, const TReal* rz, const int ldrz,
                                TReal* sz, const int ldsz)
  {
    Givens<TReal>::zero(n, rz + ldrz, ldrz, rz + ldrz + 1, ldrz ,
		           sz       , ldsz, sz        + 1, ldsz);

    while (--n > 0)
      {
        rz += ldrz + 1;
        sz += ldsz + 1;

        Givens<TReal>::zero(n, sz, ldsz, rz + ldrz + 1, ldrz,
                               sz, ldsz, sz        + 1, ldsz);
      }
  }


  template<typename TReal>
  template<typename OutputIterator>
  void
  Qrz<TReal>::bounds(int n, const TReal* const rz, const int ldrz,
                     OutputIterator out, std::vector<Givens<TReal>>& work)
  {
    const TReal* colj = rz;

    for (int j = 0; j < n; ++j, colj += ldrz)
      {
        const TReal* coli = colj + ldrz;

        for (int i = j + 1; i <= n; ++i, coli += ldrz)
          {
            TReal t = coli[j];

            for (int g = j + 1; g < i; ++g)
              {
                t = -work[g].s() * t + work[g].c() * coli[g];
              }

            work[i].gen(t, coli[i]);
          }

        *(out++) = std::abs(work[n].r());
      }
  }


  template<typename TReal>
  template<typename InputIterator>
  void
  Qrz<TReal>::permute(int n,
                      const TReal* const rz, const int ldrz, InputIterator pi,
                            TReal* const sz, const int ldsz)
  {
    for (int i = 0; i < n; ++i, ++pi)
      {
        const int j = *pi;

        std::copy_n(rz + j * ldrz, j + 1, sz + i * ldsz);
        std::fill_n(sz + i * ldsz + j + 1, n - j, 0);
      }
    std::copy_n(rz + n * ldrz, n + 1, sz + n * ldsz);

    Lapack<TReal>::geqrf(n + 1, n + 1, sz, ldsz);
  }


  template<typename TReal>
  void
  Qrz<TReal>::permute1(const int n, const int j,
                       TReal* const rz, const int ldrz)
  {
    Givens<TReal> g;

    TReal* const col0 = rz;
    TReal* const colj = rz + j * ldrz;

    for (int i = j - 1; i >= 0; --i)
      {
        g.gen(colj[i], colj[i + 1]);

        g.rot(j - i, colj - ldrz + i, -ldrz, colj - ldrz + i + 1, -ldrz);
        g.rot(n - j, colj + ldrz + i,  ldrz, colj + ldrz + i + 1,  ldrz);

        colj[i] = g.r();  colj[i + 1] = 0;
      }

    std::swap(col0[0], colj[0]);
    std::swap(col0[1], colj[1]);

    TReal* coli = col0;
    for (int i = 1; i < j; ++i)
      {
        coli += ldrz;

        g.gen(coli[i], coli[i + 1]);

        g.rot(n - i, coli + ldrz + i, ldrz, coli + ldrz + i + 1, ldrz);

        coli[i] = g.r();  coli[i + 1] = 0;
      }
  }


  template<typename TReal>
  void
  Qrz<TReal>::permute2(const int n, const int j,
                       TReal* const rz, const int ldrz)
  {
    if (j < 1)
      {
        return;
      }

    Givens<TReal> g;

    TReal* const col0 = rz;
    TReal* const colj = rz + j * ldrz;

    *(colj - ldrz + j) = 0;

    g.gen(colj[j - 1], colj[j]);

    g.rot(1    , colj - ldrz + j - 1, -ldrz, colj - ldrz + j, -ldrz,
                 colj        + j - 1, -ldrz, colj        + j, -ldrz);
    g.rot(n - j, colj + ldrz + j - 1,  ldrz, colj + ldrz + j,  ldrz);

    TReal* rii = colj - ldrz + j - 1;
    for (int i = j - 1; i > 0; --i, rii -= ldrz + 1)
      {
        *rii = 0;

        g.gen(colj[i - 1], g.r());

        g.rot(j - i + 1, colj - ldrz + i - 1, -ldrz, colj        + i, -ldrz,
                         colj        + i - 1, -ldrz, colj        + i, -ldrz);
        g.rot(n - j    , colj + ldrz + i - 1,  ldrz, colj + ldrz + i,  ldrz);
      }

    col0[0] = g.r();
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs

#endif
