/**
 * @file discontinuousgalerkin1d.cc
 * @brief NPDE homework "DiscontinuousGalerkin1D" code
 * @author Oliver Rietmann
 * @date 22.05.2019
 * @copyright Developed at ETH Zurich
 */

#include "discontinuousgalerkin1d.h"

#include <fstream>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace DiscontinuousGalerkin1D {

/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<double> compBmat(int Ml, int Mr, double h) {
  const int N = 2 * (Ml + Mr + 1);
  Eigen::SparseMatrix<double> A(N, N);
  /* SOLUTION_BEGIN */
  // Tell Eigen that the matrix is diagonal, which allows efficient
  // initialization
  A.reserve(Eigen::VectorXi::Ones(N));
  for (int i = 0; i < N; i += 2) {
    A.insert(i, i) = h;
    A.insert(i + 1, i + 1) = h * h * h / 12.0;
  }
  /* SOLUTION_END */
  return A;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double Feo(double v, double w) {
  double result;
  /* SOLUTION_BEGIN */
  auto f = [](double u) { return u * (1.0 - u); };
  auto I = [](double s) {
    double y = 0.25 - s * (1.0 - s);
    return s < 0.5 ? y : -y;
  };
  result = 0.5 * (f(v) + f(w)) - 0.5 * (I(v) - I(w));
  /* SOLUTION_END */
  return result;
}
/* SAM_LISTING_END_2 */

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

/* SAM_LISTING_BEGIN_3 */
void solveTrafficFlow() {
  int Ml = 40;
  int Mr = 40;
  int N_half = Mr + Ml + 1;
  int N = 2 * N_half;

  double h = 0.05;
  double tau = h / 3;
  double T = 1.0;
  unsigned int m = (unsigned int)(T / tau);

  /* SOLUTION_BEGIN */
  // Spatial mesh
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N_half, -2.0, 2.0);
  // Initial condition according to \prbeqref{eq:u0}. Note that we just
  // initialize the cell averages; a more refined implementation would also
  // sample the local first moments of $u_0$.
  Eigen::VectorXd mu0 = Eigen::VectorXd::Zero(N);
  auto u0 = [](double x) { return 0.0 <= x && x <= 1.0 ? 1.0 : 0.0; };
  for (int i = 0; i < N_half; ++i) {
    mu0(2 * i) = u0(x(i));
  }
  // Flux function
  auto f = [](double u) { return u * (1.0 - u); };
  // Perform fully discrete evolution
  Eigen::VectorXd mu = dgcl(mu0, f, Feo, T, Ml, Mr, h, m);
  // Retrieve cell averages 
  Eigen::VectorXd u(N_half);
  for (int i = 0; i < N_half; ++i) {
    u(i) = mu(2 * i);
  }

  std::ofstream file;
  file.open("solution.csv");
  file << x.transpose().format(CSVFormat) << std::endl;
  file << u.transpose().format(CSVFormat) << std::endl;
  file.close();
  /* SOLUTION_END */
}
/* SAM_LISTING_END_3 */

}  // namespace DiscontinuousGalerkin1D
