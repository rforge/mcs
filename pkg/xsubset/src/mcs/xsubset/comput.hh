#ifndef _MCS_XSUBSET_COMPUT_HH_
#define _MCS_XSUBSET_COMPUT_HH_


#include "xsubset.hh"


namespace MCS
{

  namespace xsubset
  {


    double
    sqr(double dx);

    double
    sign(double dx);

    double
    copysign(double dx, double dy);



    struct t_rot
    {
      double r;
      double c;
      double s;
    };


    t_rot
    rot_gen(double dx, double dy);



    struct t_chol
    {
      t_size size;
      t_size mark;
      t_dmatrix rz;
    };


    void
    chol_init(t_chol& c, t_size size);

    void
    chol_init(t_chol& c, const t_dmatrix& ay);

    void
    chol_copy(const t_chol& c, t_chol& cc);

    void
    chol_swap(t_chol& c, t_chol& cc);

    double
    chol_rss(const t_chol& c);

    void
    chol_rss(const t_chol& c, t_dbl_iter rss);

    void
    chol_shift(t_chol& c);

    void
    chol_shift(t_chol& c, t_size n);

    void
    chol_shift(const t_chol& c, t_chol& cc);

    void
    chol_drop(t_chol& c);

    void
    chol_drop(const t_chol& c, t_chol& cc);

    void
    chol_bounds(const t_chol& c, t_dbl_iter rss);

    void
    chol_permute(t_chol& c, t_size_iter index, t_dmatrix& work);



    extern double DBL_MAX;


  }

}


#endif
