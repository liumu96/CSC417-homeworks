#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>
#include <Eigen/SparseCholesky>
#include <cmath>

static const double alpha_max = 1.0;
static const double p = 0.5;
static const double c = 1e-8;
static const double cost_tol = 1e-8;
// line search is overall less costy than the full newton method, so a lower tolerance is given
static const double line_search_tol = 1e-12;

/**
 * Input:
 * @param x0 - initial point for newtons search
 * @param f(x) - function that evaluates and returns the cost function at x
 * @param g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
 * @param H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
 * @param max steps - the maximum newton iterations to take
 * @param tmp_g and tmp_H are scratch space to store gradients and hessians
 * Output:
 * @param x0 - update x0 to new value
 */
template <typename Objective, typename Jacobian, typename Hessian>
double newtons_method(
    Eigen::VectorXd &x0,
    Objective &f,
    Jacobian &g,
    Hessian &H,
    unsigned int maxSteps,
    Eigen::VectorXd &tmp_g,
    Eigen::SparseMatrixd &tmp_H)
{

   double error;
   // Begin Newton's method iteration
   for (unsigned int newton_step = 0; newton_step < maxSteps; newton_step++)
   {
      // Calculate convergence error for the current guess
      tmp_g.setZero();
      g(tmp_g, x0);
      error = tmp_g.norm();
      if (error < cost_tol)
      {
         // Enters here if current guess is good enough
         return error;
      }

      // Current guess does not converge, perhaps too ambitious
      // Proceed to line search to tune it down
      tmp_H.setZero();
      H(tmp_H, x0);
      // Quadratic minimization: H_i * d = -g_i
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
      solver.compute(tmp_H);
      Eigen::VectorXd d = solver.solve(-tmp_g);

      // Begin line search iteration
      double alpha = alpha_max;
      const double line_search_obj = f(x0) + c * d.dot(tmp_g);

      // While decrease is not sufficient and we still have space to further shrink alpha
      while (
          (f(x0 + alpha * d) > line_search_obj) && (alpha >= line_search_tol))
      {
         alpha *= p;
      }

      x0 += alpha * d;
   }

   // get error of the last step if things fall through
   error = tmp_g.norm();
   return error;
}