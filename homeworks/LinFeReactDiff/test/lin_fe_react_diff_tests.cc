#include <gtest/gtest.h>
#include <boost/filesystem.hpp>
#include "../mysolution/lin_fe_react_diff.h"

namespace LinFeReactDiff::test {

TEST(LinFeReactDiff, TestSolveFe) {
  boost::filesystem::path here = __FILE__;
  auto square_path = here.parent_path().parent_path() / "meshes/square.msh";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), square_path.string());
  auto mesh = reader.mesh();
  Eigen::VectorXd mu = solveFE(mesh);
  ASSERT_NEAR(mu(mu.size() - 1), 0.00254543, 0.00001);
  ASSERT_NEAR(mu(mu.size() - 2), 0.00145535, 0.00001);
}

TEST(LinFeReactDiff, TestEnergy) {
  boost::filesystem::path here = __FILE__;
  auto square_path = here.parent_path().parent_path() / "meshes/square.msh";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), square_path.string());
  auto mesh = reader.mesh();

  // Implementation from solution to not depend on the first task
  auto zero = [](Eigen::Vector2d x) -> double { return 0.; };
  lf::uscalfe::MeshFunctionGlobal mf_zero{zero};
  auto identity = [](Eigen::Vector2d x) -> double { return 1.; };
  lf::uscalfe::MeshFunctionGlobal mf_identity{identity};
  auto c = [](Eigen::Vector2d x) -> double { return x[0] * x[1]; };
  lf::uscalfe::MeshFunctionGlobal mf_c{c};

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::mesh::Mesh &mesh_p{*(fe_space->Mesh())};
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  const lf::base::size_type N_dofs(dofh.NoDofs());
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(mf_identity), decltype(mf_zero)>
      elmat_builder(fe_space, mf_identity, mf_zero);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_c)>
      elvec_builder(fe_space, mf_c);

  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<double>>
      rsf_edge_p = fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
  LF_ASSERT_MSG(rsf_edge_p != nullptr, "FE specification for edges missing");
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
  auto ess_bdc_flags_values_findest{
      lf::uscalfe::InitEssentialConditionFromFunction(
          dofh, *rsf_edge_p,
          [&bd_flags](const lf::mesh::Entity &edge) -> bool {
            return bd_flags(edge);
          },
          mf_zero)};

  lf::assemble::fix_flagged_solution_components<double>(
      [&ess_bdc_flags_values_findest](glb_idx_t gdof_idx) {
        return ess_bdc_flags_values_findest[gdof_idx];
      },
      A, phi);

  Eigen::VectorXd mu;
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  mu = solver.solve(phi);

  // compute energy
  double energy = computeEnergy(mesh, mu);
  ASSERT_NEAR(energy, 0.0105153, 0.00001);
}

}  // namespace LinFeReactDiff::test
