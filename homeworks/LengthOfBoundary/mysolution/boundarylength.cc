/**
 * @ file boundarylength.cc
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include <boost/filesystem/path.hpp>
#include "boundarylength.h"

namespace LengthOfBoundary {

/* SAM_LISTING_BEGIN_1 */
double volumeOfDomain(const std::shared_ptr<lf::mesh::Mesh> mesh) {
  double volume = 0.0;
  /* BEGIN_SOLUTION */
  // We loop over all entities of co-dim 0, i.e. cells and sum up their volumes (areas in 2D) via their geometry
  for (const auto& cell : mesh->Entities(0)) {
    auto cell_geo = cell.Geometry();
    volume += lf::geometry::Volume(*cell_geo);
  }
  /* END_SOLUTION */
  return volume;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double lengthOfBoundary(const std::shared_ptr<lf::mesh::Mesh> mesh) {
  double length = 0.0;
  /* BEGIN_SOLUTION */
  // Loop over all edges (co-dim 1 in 2D) and add their length if they are flagged to be on the boundary
  auto flagged_boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh, 1);
  for (const auto& entity : mesh->Entities(1)) {
    if (flagged_boundary(entity)) {
      length += lf::geometry::Volume(*entity.Geometry());
    }
  }
  /* END_SOLUTION */
  return length;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::pair<double, double> measureDomain(std::string msh_file_name) {
  double volume, length;
  /* BEGIN_SOLUTION */
  boost::filesystem::path here = __FILE__;
  auto smiley_path = here.parent_path() / msh_file_name;
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), smiley_path.string());
  auto mesh = reader.mesh();

  volume = volumeOfDomain(mesh);
  length = lengthOfBoundary(mesh);
  /* END_SOLUTION */

  return {volume, length};
}
/* SAM_LISTING_END_3 */

}  // namespace LengthOfBoundary
