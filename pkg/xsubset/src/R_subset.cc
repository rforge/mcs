/**
 * @file R_subset.cc
 *
 * @author Marc Hofmann
 */

// DEBUG
#include <iostream>
#include "mcs/core/numeric/io.hh"
// DEGUG

#include <cstring>

#include "mcs/subset/rss_selector.hh"
#include "mcs/subset/crit_selector.hh"
#include "mcs/subset/criterion.hh"



template<typename Value>
void
subset(const int observation_count, const int regressor_count, const double* const ay,
       const int mark, const double* const tolerances, const int preordering_radius,
       double* values, int* which, int& node_count)
{
  typedef typename mcs::subset::rss_selector<Value> Selector;
  typedef typename Selector::size_type Size;
  typedef typename Selector::matrix_type Matrix;
  typedef typename Selector::subset_type Subset;

  const Size m = observation_count;
  const Size n = regressor_count;
  const Size k = mark;

  const Matrix ay_mat(m, n + 1, ay);

  Selector stor(ay_mat, k);
  stor.set_preordering_radius(preordering_radius);
  stor.load_tolerances(tolerances);
  stor.select();

  which += k * n;
  for (Size j = k; j < n; ++j, which += n)
    {
      const Subset s = stor.get_best_subset(j + 1);
      values[j] = s.get_value();
      std::copy(s.begin(), s.end(), which);
    }

  node_count = stor.get_node_count();
}


template<typename Value,
	 typename Crit>
void
select(const int observation_count, const int regressor_count,
       const double* const ay, const int mark, const double tolerance,
       const int preordering_radius, double& value, int* const which,
       int& node_count, const Crit crit = Crit())
{
  typedef typename mcs::subset::crit_selector<Value, Crit> Selector;
  typedef typename Selector::size_type Size;
  typedef typename Selector::matrix_type Matrix;
  typedef typename Selector::subset_type Subset;

  const Size m = observation_count;
  const Size n = regressor_count;
  const Size k = mark;

  const Matrix ay_mat(m, n + 1, ay);

  Selector stor(ay_mat, k, crit);
  stor.set_preordering_radius(preordering_radius);
  stor.set_tolerance(tolerance);
  stor.select();

  const Subset best = stor.get_best_subset();
  value = best.get_value();
  std::copy(best.begin(), best.end(), which);

  node_count = stor.get_node_count();
}


template<typename Value>
void
select_helper(const char* criterion, const double penalty,
	      const int observation_count, const int regressor_count,
	      const double* const ay, const int mark, const double tolerance,
	      const int preordering_radius, double& value,
	      int* const which, int& node_count)
{
  typedef typename mcs::subset::aic Aic;
  typedef typename mcs::subset::aicc Aicc;
  typedef typename mcs::subset::aicc2 Aicc2;
  typedef typename mcs::subset::aicu Aicu;
  typedef typename mcs::subset::bic Bic;
  typedef typename mcs::subset::cp Cp;
  typedef typename mcs::subset::maic<double> Maic;

  if (std::strcmp(criterion, "aic") == 0)
    {
      select<Value, Aic>(observation_count, regressor_count, ay,
			 mark, tolerance, preordering_radius,
			 value, which, node_count);
    }
  else if (std::strcmp(criterion, "aicc") == 0)
    {
      select<Value, Aicc>(observation_count, regressor_count, ay,
			  mark, tolerance, preordering_radius,
			  value, which, node_count);
    }
  else if (std::strcmp(criterion, "aicc2") == 0)
    {
      select<Value, Aicc2>(observation_count, regressor_count, ay,
			   mark, tolerance, preordering_radius,
			   value, which, node_count);
    }
  else if (std::strcmp(criterion, "aicu") == 0)
    {
      select<Value, Aicu>(observation_count, regressor_count, ay,
			  mark, tolerance, preordering_radius,
			  value, which, node_count);
    }
  else if (std::strcmp(criterion, "bic") == 0)
    {
      select<Value, Bic>(observation_count, regressor_count, ay,
			 mark, tolerance, preordering_radius,
			 value, which, node_count);
    }
  else if (std::strcmp(criterion, "cp") == 0)
    {
      select<Value, Cp>(observation_count, regressor_count, ay,
			mark, tolerance, preordering_radius,
			value, which, node_count);
    }
  else if (std::strcmp(criterion, "maic") == 0)
    {
      const Maic crit(penalty);
      select<Value, Maic>(observation_count, regressor_count, ay,
			  mark, tolerance, preordering_radius,
			  value, which, node_count, crit);
    }

}


extern "C"
void
R_subset(const char** const precision, const int* const observation_count,
	 const int* const regressor_count, const double* const ay,
	 const int* const mark, const double* const tolerances,
	 const int* const preordering_radius, double* const values, int* const which,
	 int* const node_count)
{
  if (std::strcmp(*precision, "single") == 0)
    {
      subset<float>(*observation_count, *regressor_count, ay,
		    *mark, tolerances, *preordering_radius,
		    values, which, *node_count);
    }
  else if (std::strcmp(*precision, "double") == 0)
    {
      subset<double>(*observation_count, *regressor_count, ay,
		     *mark, tolerances, *preordering_radius,
		     values, which, *node_count);
    }
}


extern "C"
void
R_select(const char** const precision, const char** const criterion,
	 const double* const penalty, const int* const observation_count,
	 const int* const regressor_count, const double* const ay, const int* const mark,
	 const double* const tolerance, const int* const preordering_radius,
	 double* const value,  int* const which,  int* const node_count)
{
  if (std::strcmp(*precision, "single") == 0)
    {
      select_helper<float>(*criterion, *penalty, *observation_count,
			   *regressor_count, ay, *mark, *tolerance,
			   *preordering_radius, *value, which, *node_count);
    }
  else if (std::strcmp(*precision, "double") == 0)
    {
      select_helper<double>(*criterion, *penalty, *observation_count,
			    *regressor_count, ay, *mark, *tolerance,
			    *preordering_radius, *value, which, *node_count);
    }
}
