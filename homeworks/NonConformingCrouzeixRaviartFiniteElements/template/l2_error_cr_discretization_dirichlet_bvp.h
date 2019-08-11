/**
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss
 * @date   18.03.2019
 * @copyright Developed at ETH Zurich
 */

#ifndef NUMPDE_L2_ERROR_CR_DISCRETIZATION_DIRICHLET_BVP_H
#define NUMPDE_L2_ERROR_CR_DISCRETIZATION_DIRICHLET_BVP_H

#include <lf/io/io.h>
#include <lf/mesh/mesh.h>
#include <cmath>
#include <string>
#include "compute_cr_l2_error.h"
#include "cr_fe_space.h"
#include "solve_cr_dirichlet_bvp.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

/* SAM_LISTING_BEGIN_1 */
double L2errorCRDiscretizationDirichletBVP(const std::string &filename) {
  double l2_error;

  // TODO: task 2-14.x)
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */

  return l2_error;
}
/* SAM_LISTING_END_1 */

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // NUMPDE_L2_ERROR_CR_DISCRETIZATION_DIRICHLET_BVP_H
