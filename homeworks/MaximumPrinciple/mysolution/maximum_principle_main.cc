/**
 * @file  maximum_principle_main.cc
 * @brief NPDE homework "MaximumPrinciple" code
 * @author Oliver Rietmann
 * @date 25.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "maximum_principle.h"

#include <functional>
#include <iostream>

#include <Eigen/SparseCholesky>

using namespace MaximumPrinciple;

int main() {
  /* SAM_LISTING_BEGIN_3 */
  unsigned int M = 4;
  double c = 0.99;
  std::function<double(double, double)> f = [](double x, double y) {
    double ret = 0.0;
    /* SOLUTION_BEGIN */
    /* SOLUTION_END */
    return ret;
  };

  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  Eigen::VectorXd phi = computeLoadVector(M, f);

  Eigen::SparseMatrix<double> A = computeGalerkinMatrix(M, c);
  Eigen::VectorXd mu = solver.compute(A).solve(phi);
  std::cout << "mu = " << std::endl << mu << std::endl;
  /* SAM_LISTING_END_3 */
  Eigen::SparseMatrix<double> A_TR = computeGalerkinMatrixTR(M, c);
  Eigen::VectorXd mu_TR = solver.compute(A_TR).solve(phi);
  std::cout << "mu_TR = " << std::endl << mu_TR << std::endl;

  return 0;
}
