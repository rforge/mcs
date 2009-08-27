#ifndef _MCS_XSUBSET_XSUBSET_HH_
#define _MCS_XSUBSET_XSUBSET_HH_


#include <string>
#include <vector>
#include <functional>

#include "../xutil/vector.hh"
#include "../xutil/matrix.hh"


namespace MCS
{

  namespace xsubset
  {


    typedef std::string          t_string;
    typedef std::greater<double> t_dbl_gt;

    typedef std::size_t                  t_size;
    typedef std::vector<t_size>          t_size_array;
    typedef t_size_array::iterator       t_size_iter;
    typedef t_size_array::const_iterator t_size_iter_c;

    typedef std::vector<double>   t_dbl_array;
    typedef t_dbl_array::iterator t_dbl_iter;


    typedef xutil::vector<double> t_dvector;
    typedef xutil::matrix<double> t_dmatrix;



    struct t_args
    {
      struct
      {
        t_dmatrix  ay;
        t_size     mark;
        t_size     prad;
        t_dbl_iter tau;
      } in;

      struct
      {
        t_dbl_iter  rsel;
        t_size_iter isel;
        t_size      nvis;
      } out;
    };



    void
    select(t_dmatrix  data, t_size      mark,
           t_size     prad, t_dbl_iter  tau,
           t_dbl_iter rsel, t_size_iter isel,
           t_size&    nvis);

    void
    select(t_args& args);



  }

}


#endif
