#ifndef _MCS_XSUBSET_ENVIRON_CC_
#define _MCS_XSUBSET_ENVIRON_CC_


#include <algorithm>

#include "../xutil/algo.hh"

#include "xsubset.hh"
#include "environ.hh"
#include "comput.hh"


namespace MCS
{

  namespace xsubset
  {


    void
    node_init(t_node& n, const t_dmatrix& ay)
    {
      const t_size nsel = ay.ncol() - 1;

      n.nsel = nsel;
      n.mark = 0;
      n.isel.resize(nsel);
      for (t_size i = 0; i < nsel; ++i)
        {
          n.isel[i] = i;
        }
      chol_init(n.chol, ay);
    }


    void
    node_init(t_node& n, t_size nvar)
    {
      n.nsel = 0;
      n.mark = 0;
      n.isel.reserve(nvar - 1);
      chol_init(n.chol, nvar);
    }


    void
    node_swap(t_node& n1, t_node& n2)
    {
      std::swap(n1.nsel, n2.nsel);
      std::swap(n1.mark, n2.mark);
      n1.isel.swap(n2.isel);
      chol_swap(n1.chol, n2.chol);
    }


    void
    node_bounds(const t_node& n, t_dbl_iter bounds)
    {
      chol_bounds(n.chol, bounds);
    }


    void
    node_permute(t_node& n, t_size_iter index, t_dmatrix& work)
    {
      const t_size nsel = n.nsel;
      const t_size mark = n.mark;

      xutil::permute(n.isel.begin() + mark, nsel - mark, index);
      chol_permute(n.chol, index, work);
    }


    double
    node_rss(const t_node& n)
    {
      return chol_rss(n.chol);
    }


    void
    node_rss(const t_node& n, t_dbl_iter rss)
    {
      chol_rss(n.chol, rss);
    }


    void
    node_shift(t_node& n, t_size s)
    {
      n.mark += s;
      chol_shift(n.chol, s);
    }


    void
    node_shift(t_node& n)
    {
      ++(n.mark);
      chol_shift(n.chol);
    }


    void
    node_drop(const t_node& n1, t_node& n2)
    {
      const t_size nsel = n1.nsel;
      const t_size mark = n1.mark;

      n2.nsel = nsel - 1;
      n2.mark = mark;

      t_size_iter_c src = n1.isel.begin();
      t_size_iter dst = n2.isel.begin();
      dst = std::copy(src, src + mark, dst);
      std::copy(src + mark + 1, src + nsel, dst);

      chol_drop(n1.chol, n2.chol);
    }


    void
    env_init(t_env& env, t_args& args)
    {
      const t_size nvar = args.in.ay.ncol();
      const t_size mark = args.in.mark;

      /*
       * Root node.
       */
      t_node root;
      node_init(root, args.in.ay);
      node_shift(root, mark);

      /*
       * Preordering.
       */
      env.args.prad = (nvar - mark - 1) - args.in.prad;

      /*
       * Heuristic.
       */
      env.args.tau.resize(nvar - 1);
      for (int i = 0; i < (nvar - 1); ++i)
        {
          env.args.tau[i] = args.in.tau[i] + 1.0L;
        }

      /*
       * Regression table.
       */
      env.args.regs.resize(nvar);
      t_dbl_iter ri = args.out.rsel;
      t_size_iter si = args.out.isel;
      for (t_size i = 0; i < nvar; ++i)
        {
          *ri = DBL_MAX;

          env.args.regs[i].rsel = ri;
          env.args.regs[i].isel = si;

          ++ri;
          si += i;
        }

      /*
       * Node stack.
       */
      env.state.stack.resize(nvar);
      for (t_size i = 0; i < nvar; ++i)
        {
          node_init(env.state.stack[i], nvar);
        }
      env.state.base = env.state.stack.begin();
      env.state.top = env.state.base;
      node_swap(*((env.state.top)++), root);

      /*
       * Node count.
       */
      env.args.nvis = 0U;

      /*
       * Workspace.
       */
      node_init(env.state.node, nvar);
      env.temp.rss.resize(nvar);
      env.temp.index.resize(nvar);
      env.temp.work = t_dmatrix(nvar, nvar);
    }


    void
    env_visit(t_env& env)
    {
      ++(env.args.nvis);

      // preorder observations
      env_preord(env);

      // report regression
      env_report(env);

      // cutting test
      env_cut(env);

      while (env.state.node.mark < env.state.icut)
        {
          // compute child node
          env_push(env);

          // shift
          env_shift(env);
        }
    }


    void
    env_preord(t_env& env)
    {
      const t_size nsel = env.state.node.nsel;
      const t_size mark = env.state.node.mark;
      const t_size n = nsel - mark;

      if (n <= env.args.prad)
        {
          return;
        }

      t_dbl_iter rss = env.temp.rss.begin();
      node_bounds(env.state.node, rss);

      t_size_iter index = env.temp.index.begin();
      for (t_size i = 0; i < n; ++i)
        {
          index[i] = i;
        }
      xutil::sort_index(rss, n, index, t_dbl_gt());

      node_permute(env.state.node, index, env.temp.work);
    }


    void
    env_report(t_env& env)
    {
      const t_size nsel = env.state.node.nsel;
      const t_size mark = env.state.node.mark;
      const t_size n = nsel - mark;

      t_dbl_iter rss = env.temp.rss.begin();
      node_rss(env.state.node, rss);

      t_reg_iter regs = env.args.regs.begin();
      t_size_iter isel = env.state.node.isel.begin();
      for (t_size i = 0; i < n; ++i)
        {
          if (rss[i] < *(regs[mark + i + 1].rsel))
            {
              *(regs[mark + i + 1].rsel) = rss[i];
              xutil::copy_n(isel, mark + i + 1, regs[mark + i + 1].isel);
            }
        }
    }


    void
    env_cut(t_env& env)
    {
      const t_size nsel = env.state.node.nsel;
      const t_size mark = env.state.node.mark;

      t_dbl_iter tau = env.args.tau.begin();
      t_reg_iter regs = env.args.regs.begin();

      int i;
      for (i = nsel - 1; i > mark; --i)
        {
          double b = tau[i] * env.state.rss;
          if (b < *(regs[i].rsel))
            {
              break;
            }
        }
      env.state.icut = i;
    }


    void
    env_shift(t_env& env)
    {
      node_shift(env.state.node);
    }


    void
    env_push(t_env& env)
    {
      node_drop(env.state.node, *((env.state.top)++));
    }


    void
    env_pop(t_env& env)
    {
      node_swap(env.state.node, *(--(env.state.top)));
      env.state.rss = node_rss(env.state.node);
    }


  }

}


#endif
