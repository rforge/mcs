#ifndef _MCS_XSUBSET_ENVIRON_HH_
#define _MCS_XSUBSET_ENVIRON_HH_


#include <vector>
#include <list>

#include "xsubset.hh"
#include "comput.hh"


namespace MCS
{

  namespace xsubset
  {


    struct t_reg
    {
      t_dbl_iter  rsel;
      t_size_iter isel;
    };

    typedef std::vector<t_reg>    t_reg_array;
    typedef t_reg_array::iterator t_reg_iter;



    struct t_node
    {
      t_size       nsel;
      t_size       mark;
      t_size_array isel;
      t_chol       chol;
    };

    typedef std::vector<t_node>    t_node_array;
    typedef t_node_array::iterator t_node_iter;


    void
    node_init(t_node& n,
              const t_dmatrix& ay);

    void
    node_init(t_node& n,
              t_size nsel);

    void
    node_swap(t_node& n1,
              t_node& n2);

    void
    node_bounds(const t_node& n,
                t_dbl_iter bounds);

    void
    node_permute(t_node& n,
                 t_size_iter index,
                 t_dmatrix& work);

    double
    node_rss(const t_node& n);

    void
    node_rss(const t_node& n,
             t_dbl_iter rss);

    void
    node_report(const t_node& n,
                t_reg_iter regs);

    void
    node_shift(t_node& n,
               t_size s);

    void
    node_shift(t_node& n);

    void
    node_drop(const t_node& n1,
              t_node& n2);



    struct t_env
    {
      struct
      {
        t_size       prad;
        t_dbl_array  tau;
        t_reg_array  regs;
        unsigned int nvis;
      } args;

      struct
      {
        t_node_array stack;
        t_node_iter  base;
        t_node_iter  top;

        t_node node;
        double rss;
        t_size icut;
      } state;

      struct
      {
        t_dbl_array  rss;
        t_size_array index;
        t_dmatrix    work;
      } temp;
    };


    void
    env_init(t_env& env,
             t_args& args);

    void
    env_visit(t_env& env);

    void
    env_preord(t_env& env);

    void
    env_report(t_env& env);

    void
    env_cut(t_env& env);

    void
    env_shift(t_env& env);

    void
    env_push(t_env& env);

    void
    env_pop(t_env& env);


  }

}


#endif
