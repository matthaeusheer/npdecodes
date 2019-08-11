/**
 This homework problem consists of reading a simple, gmesh generated, mesh on
 the unit square and solving a simple reaction diffusion system using LehrFEM++
 */

#include "lin_fe_react_diff.h"

namespace LinFeReactDiff {

/**
 * @brief generate hierarchy of meshes
 * @param levels: number of refinement steps
 */

std::shared_ptr<lf::refinement::MeshHierarchy> generateMeshHierarchy(
    const lf::base::size_type levels) {
  // set path
  boost::filesystem::path here = __FILE__;
  auto square_path = here.parent_path().parent_path() / "meshes/square.msh";

  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), square_path.string());
  auto mesh = reader.mesh();

  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh, levels);

  return multi_mesh_p;
}

/**
 * @brief solve the linear variational problem
 *        int_{\Omgea} grad(u)grad(v) dx = \int_{\Omgea} cv dx for all v
 *        with 0 Dirichlet boundary condition.
 * @param mesh: mesh discretization of computational Domain \Omega
 */
Eigen::VectorXd solveFE(std::shared_ptr<const lf::mesh::Mesh> mesh) {
  // Initialize mesh functions for solving the BVP:
  // \int_{\Omega} grad(u)grad(v) + uv dx = \_int(\Omega}cv dx
  // with Dirichlet boundary conditions fixed to 0.
  // used for the boundary conditions
  auto zero = [](Eigen::Vector2d x) -> double { return 0.; };
  lf::uscalfe::MeshFunctionGlobal mf_zero{zero};

  // used for the Galerkin Matrix
  auto identity = [](Eigen::Vector2d x) -> double { return 1.; };
  lf::uscalfe::MeshFunctionGlobal mf_identity{identity};

  // used for the load vector
  auto c = [](Eigen::Vector2d x) -> double { return x[0] * x[1]; };
  lf::uscalfe::MeshFunctionGlobal mf_c{c};

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::mesh::Mesh &mesh_p{*(fe_space->Mesh())};

  // Initialize dof handler
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  const lf::base::size_type N_dofs(dofh.NoDofs());

  // Set up Galerkin Matrix
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);

  // Populate Galerkin Matrix
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(mf_identity), decltype(mf_zero)>
      elmat_builder(fe_space, mf_identity, mf_zero);

  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // Set up load vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_c)>
      elvec_builder(fe_space, mf_c);

  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<double>>
      rsf_edge_p = fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
  LF_ASSERT_MSG(rsf_edge_p != nullptr, "FE specification for edges missing");

  // Fetch flags and values for degrees of freedom located on Dirichlet
  // edges.
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
  auto ess_bdc_flags_values_findest{
      lf::uscalfe::InitEssentialConditionFromFunction(
          dofh, *rsf_edge_p,
          [&bd_flags](const lf::mesh::Entity &edge) -> bool {
            return bd_flags(edge);
          },
          mf_zero)};

  // Eliminate Dirichlet dofs from linear system
  lf::assemble::fix_flagged_solution_components<double>(
      [&ess_bdc_flags_values_findest](glb_idx_t gdof_idx) {
        return ess_bdc_flags_values_findest[gdof_idx];
      },
      A, phi);

  // Solve System
  Eigen::VectorXd mu;
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return mu;
}

double computeEnergy(std::shared_ptr<const lf::mesh::Mesh> mesh,
                     Eigen::VectorXd mu) {
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);

  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  const lf::base::size_type N_dofs(dofh.NoDofs());

  auto identity = [](Eigen::Vector2d x) -> double { return 1.; };
  lf::uscalfe::MeshFunctionGlobal mf_identity{identity};

  auto zero = [](Eigen::Vector2d x) -> double { return 0.; };
  lf::uscalfe::MeshFunctionGlobal mf_zero{zero};

  // Matrix in triplet format holding Stiffness matrix.
  lf::assemble::COOMatrix<double> Stiffness(N_dofs, N_dofs);
  // Assemble Stiffness matrix for \int_{\Omega} uv dx
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */

  Eigen::SparseMatrix<double> Stiffness_crs = Stiffness.makeSparse();

  // Matrix in triplet format holding Mass matrix.
  lf::assemble::COOMatrix<double> Mass(N_dofs, N_dofs);
  // Assemble Mass matrix for \int_{\Omega} uv dx
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  double energy_stiffness_sq;
  double energy_mass_sq;
  // energy_stiffness_sq = 1' A 1
  // energy_mass_sq = \mu' M \mu

  /* SOLUTION_BEGIN */
  /* SOLUTION_END */

  return std::sqrt(energy_stiffness_sq + energy_mass_sq);
}

}  // namespace LinFeReactDiff
