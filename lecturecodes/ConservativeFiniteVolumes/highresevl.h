#ifndef HIGHRESEVL_HPP
#define HIGHRESEVL_HPP

#include <Eigen/Dense>

#include "ode.h"
#include "slopelimfluxdiff.h"

// \ref{cpp:highresevl} in lecture document
/* SAM_LISTING_BEGIN_1 */
// Arguments:
//   real numbers \Blue{$a,b$}. \Blue{$a<b$}, the boundaries of the domain
//   unsigned int \Blue{$N$} the number of grid cells
//   Functor \Blue{$u0 : \mathbb R \mapsto \mathbb R$}: initial data
//   real number \Blue{$T>0$}: final time
//   Functor \Blue{$F=F(v,w)$}:  2-point numerical flux function
//   Functor \Blue{$slopes = \sigma(v,u,w)$}: 3-point slope recostruction rule
//  (Note: no division by \Blue{$h$} needs to be done in slopt computation
//
// returns:
//   Vector of cell averages at final time
//
// finite volume discrete evolution in conservation form with linear
// reconstruction, see \eqref{eq:lccf}.
// Cauchy problem over time \Blue{$[0,T]$}, equidistant mesh
template <typename FunctionU0, typename FunctionF, typename FnSlopes>
Eigen::VectorXd highresevl(double a, double b, unsigned N, FunctionU0 u0,
                           double T, FunctionF F, FnSlopes slopes) {
  double h = (b - a) / N;  // mesh width
  // positions of grid points
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N, a + 0.5 * h, b - 0.5 * h);
  // vector of initial cell averages (column vector) from sampling \Blue{$u_0$}
  Eigen::VectorXd mu0 = h * x.unaryExpr(u0);

  // right hand side lambda function for ODE solver
  auto odefun = [&](const Eigen::VectorXd& mu, Eigen::VectorXd& dmdt,
                    double t) {
    dmdt = -1. / h * slopelimfluxdiff<FunctionF, FnSlopes>(h, mu, F, slopes);
  };

  // timestepping by explicit adaptive Runge-Kutta single-step
  // method of order 5. Adaptivity control according to \ncseref{sec:ssctrl}
  double abstol = 1E-8, reltol = 1E-6;
  Eigen::VectorXd t;   // Temporal adaptive integration mesh
  Eigen::MatrixXd MU;  // State vectors \Blue{$\vec{\mubf}^{(k)}$}
  std::tie(t, MU) = ode45(odefun, 0, T, mu0, abstol, reltol);
  // Extract state vector for final time
  return MU.col(t.size() - 1);
}
/* SAM_LISTING_END_1 */
}  // namespace ConsFV

#endif
