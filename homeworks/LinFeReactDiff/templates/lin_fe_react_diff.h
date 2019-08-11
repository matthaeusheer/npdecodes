/**
 This homework problem consists of reading a simple, gmesh generated, mesh on
 the unit square and solving a simple reaction diffusion system using LehrFEM++
 */

#include <fstream>
#include <iomanip>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <boost/filesystem.hpp>

namespace LinFeReactDiff {

using glb_idx_t = lf::assemble::glb_idx_t;

std::shared_ptr<lf::refinement::MeshHierarchy> generateMeshHierarchy(
    const lf::base::size_type levels);

Eigen::VectorXd solveFE(std::shared_ptr<const lf::mesh::Mesh> mesh);

double computeEnergy(std::shared_ptr<const lf::mesh::Mesh> mesh,
                     Eigen::VectorXd mu);

}  // namespace LinFeReactDiff
