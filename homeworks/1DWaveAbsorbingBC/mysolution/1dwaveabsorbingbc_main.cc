/**
 * @file 1dwaveabsorbingbc_main.cc
 * @brief NPDE homework "1DWaveAbsorbingBC" code
 * @author Oliver Rietmann
 * @date 08.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "1dwaveabsorbingbc.h"

#include <fstream>

#include <Eigen/Core>

using namespace WaveAbsorbingBC1D;

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  double c = 1.0;
  double T = 7.0;
  unsigned int N = 100;
  unsigned int m = 2000;
  Eigen::MatrixXd R = waveLeapfrogABC(c, T, N, m);

  std::pair<Eigen::VectorXd, Eigen::VectorXd> energies =
      computeEnergies(R, c, T / m);
  Eigen::VectorXd E_pot = energies.first;
  Eigen::VectorXd E_kin = energies.second;
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(m + 1, 0.0, T);

  // print the data, e.g. to a .csv file, in a suitable way
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */

  return 0;
}
