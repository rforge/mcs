#ifndef _MCS_XSUBSET_XSUBSET_CC_
#define _MCS_XSUBSET_XSUBSET_CC_


#include "xsubset.hh"
#include "environ.hh"


namespace MCS
{

  namespace xsubset
  {


    void
    select(t_dmatrix  ay,   t_size      mark,
           t_size     prad, t_dbl_iter  tau,
           t_dbl_iter rsel, t_size_iter isel,
           t_size&    nvis)
    {
      t_args args;

      args.in.ay = ay;
      args.in.mark = mark;
      args.in.prad = prad;
      args.in.tau = tau;
      args.out.rsel = rsel;
      args.out.isel = isel;
      args.out.nvis = 0;

      select(args);

      nvis = args.out.nvis;
    }


    void select(t_args& args)
    {
      t_env env;

      // init
      env_init(env, args);

      while (env.state.top != env.state.base)
        {
          // pop node
          env_pop(env);

          // visit current node
          env_visit(env);
        }

      args.out.nvis = env.args.nvis;
    }


  }

}


#endif
