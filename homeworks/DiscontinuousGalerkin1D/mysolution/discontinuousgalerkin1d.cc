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

Eigen::SparseMatrix<double> compBmat(int Ml, int Mr, double h) {
  const int N = 2 * (Ml + Mr + 1);
  Eigen::SparseMatrix<double> A(N, N);
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return A;
}

double Feo(double v, double w) {
  double result;
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return result;
}

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

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
  /* SOLUTION_END */
}

}  // namespace DiscontinuousGalerkin1D
