/**
 * @file
 * @brief NPDE homework ErrorEstimatesForTraces
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "tee_lapl_robin_assembly.h"

namespace ErrorEstimatesForTraces {

Eigen::VectorXd solveBVP(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space) {
  // I : Creating coefficients as Lehrfem++ mesh functions
  // Coefficients used in the class template
  // ReactionDiffusionElementMatrixProvider<SCALAR,DIFF_COEFF,REACTION_COEFF>
  auto alpha =
      lf::uscalfe::MeshFunctionGlobal([](coord_t x) -> double { return 1.0; });
  auto gamma =
      lf::uscalfe::MeshFunctionGlobal([](coord_t x) -> double { return 0.0; });
  // Coefficients used in the class template
  // MassEdgeMatrixProvider< SCALAR, COEFF, EDGESELECTOR >
  auto eta =
      lf::uscalfe::MeshFunctionGlobal([](coord_t x) -> double { return 1.0; });
  // Right-hand side source function f
  auto f = lf::uscalfe::MeshFunctionGlobal(
      [](coord_t x) -> double { return std::cos(x.norm()); });

  // pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // II: Instantiating finite element matrix for the Robin bilinear form
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NoDofs());
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Right hand side vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();  // has to be zero initially

  // III : Computing mass (volume) matrix
  // Initialize object taking care of local mass (volume) computations.
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                                      decltype(gamma)>
      elmat_builder(fe_space, alpha, gamma);
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix A.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // IV : Computing contribution from the boundary bilinear form associated to
  // the Robin boundary conditions.
  // Obtain an array of boolean flags for the edges of the mesh, 'true'
  // indicates that the edge lies on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Creating a predicate that will guarantee that the computations are carried
  // only on the edges of the mesh using the boundary flags
  auto edges_predicate = [&bd_flags](const lf::mesh::Entity& edge) -> bool {
    return bd_flags(edge);
  };
  lf::uscalfe::MassEdgeMatrixProvider<double, decltype(eta),
                                      decltype(edges_predicate)>
      edgemat_builder(fe_space, eta, edges_predicate);
  // Invoke assembly on edges by specifying co-dimension = 1
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, edgemat_builder, A);

  // V: Computing right-hand side vector
  // Assemble volume part of right-hand side vector depending on the source
  // function f. Initialize object taking care of local computations on all
  // cells.
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(f)>
      elvec_builder(fe_space, f);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // Assembly completed: Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Solve linear system using Eigen's sparse direct elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  Eigen::VectorXd sol_vec = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

  return sol_vec;
}

double bdFunctionalEval(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space,
    Eigen::VectorXd& coeff_vec) {
  double bd_functional_val = 0;

  /* BEGIN_SOLUTION */

  /* END_SOLUTION */
  return bd_functional_val;
}

}  // namespace ErrorEstimatesForTraces
