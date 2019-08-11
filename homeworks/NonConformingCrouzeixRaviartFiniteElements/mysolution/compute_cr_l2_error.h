/**
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss
 * @date   18.03.2019
 * @copyright Developed at ETH Zurich
 */

#ifndef NUMPDE_COMPUTE_CR_L2_ERROR_H
#define NUMPDE_COMPUTE_CR_L2_ERROR_H

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include "cr_fe_space.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTION>
double computeCRL2Error(std::shared_ptr<CRFeSpace> fe_space,
                        const Eigen::VectorXd &mu, FUNCTION &&u) {
  double l2_error = 0.;

  // TODO task 2-14.w)
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */

  return std::sqrt(l2_error);
}
/* SAM_LISTING_END_1 */

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // NUMPDE_COMPUTE_CR_L2_ERROR_H
