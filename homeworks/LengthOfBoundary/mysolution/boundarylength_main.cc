/**
 * @ file boundarylength_main.cc
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "boundarylength.h"

int main(int argc, char *argv[]) {

  std::cout << "Loading mesh \"square.msh\"..." << std::endl;
  auto domain_measurement = LengthOfBoundary::measureDomain("square.msh");
  std::cout << "Volume: " << domain_measurement.first << std::endl;
  std::cout << "Boundary length: " << domain_measurement.second << std::endl;
}
