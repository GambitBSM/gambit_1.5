// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef FIXED_POINT_ITERATOR_H
#define FIXED_POINT_ITERATOR_H

#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <Eigen/Core>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "logger.hpp"
#include "wrappers.hpp"
#include "error.hpp"
#include "ewsb_solver.hpp"
#include "gsl_utils.hpp"
#include "gsl_vector.hpp"

namespace flexiblesusy {

namespace fixed_point_iterator {

class Convergence_tester_absolute {
public:
   explicit Convergence_tester_absolute(double precision_ = 1.0e-2)
      : precision(precision_)
   {}

   std::string name() const { return "Convergence_tester_absolute"; }

   /**
    * Test whether the absolute value of the residual, defined by
    * \f$|a-b| = \sqrt{\sum_i (a_i - b_i)^2}\f$,
    * is less than the set precision.
    *
    * @param a GSL vector
    * @param b GSL vector
    * @return GSL error code (GSL_SUCCESS or GSL_CONTINUE)
    */
   int operator()(const GSL_vector& a, const GSL_vector& b) const {
      if (a.size() != b.size()) throw SetupError("Error: vectors have different size.");

      const auto dimension = a.size();
      double residual = 0.;

      if (precision < 0.)
         GSL_ERROR("absolute tolerance is negative", GSL_EBADTOL);

      for (std::size_t i = 0; i < dimension; ++i)
         residual += Sqr(a[i] - b[i]);

      residual = Sqrt(residual);

      return (residual < precision ? GSL_SUCCESS : GSL_CONTINUE);
   }

private:
   double precision;                 ///< precision goal
};

class Convergence_tester_relative {
public:
   explicit Convergence_tester_relative(double precision_ = 1.0e-2)
      : precision(precision_)
   {}

   std::string name() const { return "Convergence_tester_relative"; }

   /**
    * Test whether the relative difference is less than the set
    * precision. The relative difference test used here is carried out
    * by applying \a MaxRelDiff to each element of the vector.
    *
    * @param a GSL vector
    * @param b GSL vector
    * @return GSL error code (GSL_SUCCESS or GSL_CONTINUE)
    */
   int operator()(const GSL_vector& a, const GSL_vector& b) const {
      if (a.size() != b.size()) throw SetupError("Error: vectors have different size.");

      const auto dimension = a.size();
      double rel_diff = 0.;

      if (precision < 0.)
         GSL_ERROR("relative tolerance is negative", GSL_EBADTOL);

      for (std::size_t i = 0; i < dimension; ++i) {
         rel_diff = MaxRelDiff(a[i], b[i]);

         if (rel_diff > precision)
            return GSL_CONTINUE;
      }

      return GSL_SUCCESS;
   }

private:
   double precision;                 ///< precision goal
};

template <std::size_t dimension>
class Convergence_tester_tadpole {
public:
   using Vector_t = Eigen::Matrix<double,dimension,1>;
   using Function_t = std::function<Vector_t(const Vector_t&)>;

   Convergence_tester_tadpole(double precision_,
                              const Function_t& tadpole_function_)
      : precision(precision_)
      , tadpole_function(tadpole_function_)
   {}

   std::string name() const { return "Convergence_tester_tadpole"; }

   /**
    * Test whether the relative difference is less than the set
    * precision. The relative difference test used here is carried out
    * by applying \a MaxRelDiff to each element of the vector. If the
    * relative difference is below the precision, it is tested whether
    * the tadpoles are below the precision. If the tadpoles are larger
    * than the precision, GSL_CONTINUE is returned.
    *
    * @param a GSL vector
    * @param b GSL vector
    * @return GSL error code (GSL_SUCCESS or GSL_CONTINUE)
    */
   int operator()(const GSL_vector& a, const GSL_vector& b) const {
      if (a.size() != b.size()) throw SetupError("Error: vectors have different size.");

      if (precision < 0.)
         GSL_ERROR("relative tolerance is negative", GSL_EBADTOL);

      const double max_rel_diff =
         MaxRelDiff(to_eigen_vector(a), to_eigen_vector(b));

      if (max_rel_diff > precision)
         return GSL_CONTINUE;

      static const double eps = 10*std::pow(10., -std::numeric_limits<double>::digits10);

      if (max_rel_diff < eps)
         return GSL_SUCCESS;

      return check_tadpoles(a);
   }

private:
   double precision;                 ///< precision goal
   const Function_t tadpole_function; ///< function to calculate tadpole

   int check_tadpoles(const GSL_vector& x) const {
      const GSL_vector t(to_GSL_vector(tadpole_function(to_eigen_vector(x))));
      return gsl_multiroot_test_residual(t.raw(), precision);
   }
};

} // namespace fixed_point_iterator

/**
 * @class Fixed_point_iterator
 * @brief Does fixed point iteration
 * @author Dylan Harries, Alexander Voigt
 * @tparam dimension dimension of function
 * @tparam Convergence_tester function for relative comparison
 *    of subsequent iteration steps
 *
 * The user has to provide the function (of which a fixed point should
 * be found) of the type \a Function_t. This function gets as
 * arguments a Eigen vector of length \a dimension and returns a
 * vector with the next point.
 *
 * @note The standard relative convergence criterion
 * \f$\text{MaxRelDiff}(x_{n+1}, x_{n}) < \text{precision}\f$ is not
 * very good: The iteration might converge slowly.  This means, that
 * subsequent steps are very close to each other, but \f$x_n\f$ might
 * not be close to the true fixed point.
 *
 * @todo implement check for no progress towards solution
 */
template <std::size_t dimension, class Convergence_tester = fixed_point_iterator::Convergence_tester_relative>
class Fixed_point_iterator : public EWSB_solver {
public:
   using Vector_t = Eigen::Matrix<double,dimension,1>;
   using Function_t = std::function<Vector_t(const Vector_t&)>;

   Fixed_point_iterator() = default;
   template <typename F>
   Fixed_point_iterator(F&&, std::size_t, const Convergence_tester&);
   virtual ~Fixed_point_iterator() = default;

   template <typename F>
   void set_function(F&& f) { function = std::forward<F>(f); }
   void set_max_iterations(std::size_t n) { max_iterations = n; }
   int find_fixed_point(const Eigen::VectorXd&);

   // EWSB_solver interface methods
   virtual std::string name() const override { return "Fixed_point_iterator<" + convergence_tester.name() + ">"; }
   virtual int solve(const Eigen::VectorXd&) override;
   virtual Eigen::VectorXd get_solution() const override;

private:
   std::size_t max_iterations{100};         ///< maximum number of iterations
   GSL_vector xn{dimension};                ///< current iteration point
   GSL_vector fixed_point{dimension};       ///< vector of fixed point estimate
   Function_t function{nullptr};            ///< function defining fixed point
   Convergence_tester convergence_tester{}; ///< convergence tester

   int fixed_point_iterator_iterate();
   void print_state(std::size_t) const;
   static int gsl_function(const gsl_vector*, void*, gsl_vector*);
};

/**
 * Constructor
 *
 * @param function_ pointer to the function to find fixed point for
 * @param max_iterations_ maximum number of iterations
 * @param convergence_tester_ convergence tester
 */
template <std::size_t dimension, class Convergence_tester>
template <typename F>
Fixed_point_iterator<dimension,Convergence_tester>::Fixed_point_iterator(
   F&& function_,
   std::size_t max_iterations_,
   const Convergence_tester& convergence_tester_
)
   : max_iterations(max_iterations_)
   , function(std::forward<F>(function_))
   , convergence_tester(convergence_tester_)
{
}

/**
 * Start the iteration
 *
 * @param start starting point
 *
 * @return GSL error code (GSL_SUCCESS if fixed point found)
 */
template <std::size_t dimension, class Convergence_tester>
int Fixed_point_iterator<dimension,Convergence_tester>::find_fixed_point(
   const Eigen::VectorXd& start
)
{
   if (!function)
      throw SetupError("Fixed_point_iterator: function not callable");

   int status;
   std::size_t iter = 0;

#ifndef ENABLE_DEBUG
   gsl_set_error_handler_off();
#endif

   fixed_point = xn = to_GSL_vector(start);

#ifdef ENABLE_VERBOSE
   print_state(iter);
#endif

   do {
      iter++;
      status = fixed_point_iterator_iterate();

#ifdef ENABLE_VERBOSE
      print_state(iter);
#endif

      if (status)   // check if iterator has problems
         break;

      status = convergence_tester(fixed_point, xn);

   } while (status == GSL_CONTINUE && iter < max_iterations);

   VERBOSE_MSG("\t\t\tFixed_point_iterator status = "
               << gsl_strerror(status));

   return status;
}

/**
 * Perform a single step of the fixed point iteration
 *
 * @return GSL error code
 */
template <std::size_t dimension, class Convergence_tester>
int Fixed_point_iterator<dimension,Convergence_tester>::fixed_point_iterator_iterate()
{
   xn = fixed_point;

   void* parameters = &function;

   int status = gsl_function(xn.raw(), parameters, fixed_point.raw());

   if (status != GSL_SUCCESS)
      return GSL_EBADFUNC;

   // For safety, include a check for nans or infs here (which
   // should be sufficient for now)
   if (!is_finite(fixed_point))
      GSL_ERROR("update point is not finite", GSL_EBADFUNC);

   return GSL_SUCCESS;
}

/**
 * Print state of the fixed point iterator
 *
 * @param iteration iteration number
 */
template <std::size_t dimension, class Convergence_tester>
void Fixed_point_iterator<dimension,Convergence_tester>::print_state(std::size_t iteration) const
{
   VERBOSE_MSG("\t\t\tIteration n = " << iteration
               << ": x_{n} = " << xn
               << ", x_{n+1} = " << fixed_point);
}

template <std::size_t dimension, class Convergence_tester>
int Fixed_point_iterator<dimension,Convergence_tester>::gsl_function(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (!is_finite(x)) {
      gsl_vector_set_all(f, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   Function_t* fun = static_cast<Function_t*>(params);
   int status = GSL_SUCCESS;
   const Vector_t arg(to_eigen_vector(x));
   auto result = arg;

   try {
      result = (*fun)(arg);
      status = GSL_SUCCESS;
   } catch (const flexiblesusy::Error&) {
      status = GSL_EDOM;
   }

   copy(result, f);

   return status;
}

template <std::size_t dimension, class Convergence_tester>
int Fixed_point_iterator<dimension,Convergence_tester>::solve(const Eigen::VectorXd& start)
{
   return (find_fixed_point(start) == GSL_SUCCESS ?
           EWSB_solver::SUCCESS : EWSB_solver::FAIL);
}

template <std::size_t dimension, class Convergence_tester>
Eigen::VectorXd Fixed_point_iterator<dimension,Convergence_tester>::get_solution() const
{
   return to_eigen_vector(fixed_point);
}

} // namespace flexiblesusy

#endif
