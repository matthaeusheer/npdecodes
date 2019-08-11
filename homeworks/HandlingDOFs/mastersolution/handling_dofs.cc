/**
 * @file
 * @brief NPDE homework "Handling degrees of freedom (DOFs) in LehrFEM++"
 * @author Julien Gacon
 * @date March 1st, 2019
 * @copyright Developed at ETH Zurich
 */

#include "handling_dofs.h"

namespace HandlingDOFs {

/* SAM_LISTING_BEGIN_1 */
std::array<std::size_t, 3> countEntityDofs(
    const lf::assemble::DofHandler& dofhandler) {
  std::array<std::size_t, 3> entityDofs;
  /* BEGIN_SOLUTION */
  // Idea: iterate over entities in the mesh and get interior number of dofs for
  // each
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  for (std::size_t codim = 0; codim <= 2; ++codim) {
    entityDofs[codim] = 0;
    for (const auto& el : mesh->Entities(codim)) {
      if (el.RefEl() == lf::base::RefEl::kQuad()) {
        throw "Only triangular meshes are allowed!";
      }
      entityDofs[codim] += dofhandler.NoInteriorDofs(el);
    }
  }
  /* END_SOLUTION */
  return entityDofs;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::size_t countBoundaryDofs(const lf::assemble::DofHandler& dofhandler) {
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  // given an entity, bd\_flags(entity) == true, if the entity is on the
  // boundary
  lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags(
      lf::mesh::utils::flagEntitiesOnBoundary(mesh));
  std::size_t no_dofs_on_bd = 0;
  /* BEGIN_SOLUTION */
  // Edges and nodes can be on the boundary
  for (const auto& edge : mesh->Entities(1)) {
    if (bd_flags(edge)) {
      no_dofs_on_bd += dofhandler.NoInteriorDofs(edge);
    }
  }
  for (const auto& node : mesh->Entities(2)) {
    if (bd_flags(node)) {
      no_dofs_on_bd += dofhandler.NoInteriorDofs(node);
    }
  }
  /* END_SOLUTION */
  return no_dofs_on_bd;
}
/* SAM_LISTING_END_2 */

// clang-format off
/* SAM_LISTING_BEGIN_3 */
double integrateLinearFEFunction(
    const lf::assemble::DofHandler& dofhandler,
    const Eigen::VectorXd& mu) {
  double I = 0;
  /* BEGIN_SOLUTION */
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  for (const auto& cell : mesh->Entities(0)) {
    // check if we the FE space is really $\Cs_1^0$
    if (dofhandler.NoLocalDofs(cell) != 3) {
      throw "Not a S_1^0 FE space!";
    }
    // iterate over dofs
    auto int_dofs = dofhandler.GlobalDofIndices(cell);
    for (auto dof_idx_p = int_dofs.begin(); dof_idx_p < int_dofs.end();
         ++dof_idx_p) {
      // local integral of the basis function associated with this dof:
      // in linear Lagrangian FE, the integral over the basis functions over
      // a triangle K is: 1/3*vol(K)
      const double I_bary = 1.0 / 3.0 * lf::geometry::Volume(*cell.Geometry());
      // multiply by the value at the dof to get local contribution
      I += I_bary * mu(*dof_idx_p);
    }
  }
  /* END_SOLUTION */
  return I;
}
/* SAM_LISTING_END_3 */
// clang-format on

/* SAM_LISTING_BEGIN_4 */
double integrateQuadraticFEFunction(const lf::assemble::DofHandler& dofhandler,
                                    const Eigen::VectorXd& mu) {
  double I = 0;
  /* BEGIN_SOLUTION */
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  for (const auto& cell : mesh->Entities(0)) {
    // check if we the FE space is really $\Cs_2^0$
    if (dofhandler.NoLocalDofs(cell) != 6) {
      throw "Not a S_2^0 FE space!";
    }
    const double weight = 1.0 / 3.0 * lf::geometry::Volume(*cell.Geometry());
    // iterate over dofs
    auto int_dofs = dofhandler.GlobalDofIndices(cell);
    // The integrated basis functions associated with the nodes are 0:
    //  $\int_K b_K^j dx = 0$ for $j = 1,2,3$. Skip!
    for (int l = 3; l < 6; ++l) {
      // The integrated basis functions associated with the edges are:
      // $\int_K b_K^j dx = |K|/3$ for $j = 4,5,6$
      // multiply by the value at the dof to get local contribution
      I += (weight * mu(int_dofs[l]));
    }
  }
  /* END_SOLUTION */
  return I;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd convertDOFsLinearQuadratic(
    const lf::assemble::DofHandler& dofh_Linear_FE,
    const lf::assemble::DofHandler& dofh_Quadratic_FE,
    const Eigen::VectorXd& mu) {
  if (dofh_Linear_FE.Mesh() != dofh_Quadratic_FE.Mesh()) {
    throw "Underlying meshes must be the same for both DOF handlers!";
  }
  std::shared_ptr<const lf::mesh::Mesh> mesh =
      dofh_Linear_FE.Mesh();                         // get the mesh
  Eigen::VectorXd zeta(dofh_Quadratic_FE.NoDofs());  // initialise empty zeta
  // safety guard: always set zero if you're not sure to set every entry later
  // on for us this shouldn't be a problem, but just to be sure
  zeta.setZero();

  for (const auto& cell : mesh->Entities(0)) {
    // check if the spaces are actually linear and quadratic
    /* BEGIN_SOLUTION */
    if (dofh_Linear_FE.NoLocalDofs(cell) != 3 ||
        dofh_Quadratic_FE.NoLocalDofs(cell) != 6) {
      throw "dofh_Linear_FE must have 3 dofs per cell and dofh_Quadratic_FE 6!";
    }
    /* END_SOLUTION */
    // get the global dof indices of the linear and quadratic FE spaces, note
    // that the vectors obey the LehrFEM++ numbering, which we will make use of
    // lin\_dofs will have size 3 for the 3 dofs on the nodes and
    // quad\_dofs will have size 6, the first 3 entries being the nodes and
    // the last 3 the edges
    /* BEGIN_SOLUTION */
    auto lin_dofs = dofh_Linear_FE.GlobalDofIndices(cell);
    auto quad_dofs = dofh_Quadratic_FE.GlobalDofIndices(cell);
    /* END_SOLUTION */
    // assign the coefficients of mu to the correct entries of zeta, use
    // the previous subproblem 2-9.a
    /* BEGIN_SOLUTION */
    for (std::size_t l = 0; l <= 2; ++l) {
      // Let $p_1,p_2,p_3$ be the nodes of the triangle and $p_4,p_5,p_6$ the
      // edge midpoints. Let $\lambda_1,\lambda_2,\lambda_3$ be the local basis
      // functions of $\Cs_1^0$ and $b_1,..,b_6$ for $\Cs_2^0$. The values of
      // $\lambda_i$ in the interpolation nodes $p_1,..,p_6$ are the
      // coefficients of writing $\lambda_i$ as lin. comb of $b_1,..b_6$.
      // From 2-9.a we know $\lambda_l(p_l) = 1$; $l = 1,2,3$.
      // Hence we simply copy the coefficient of $\lambda_l$ to $b_l$ for
      // $l=1,2,3$.
      zeta(quad_dofs[l]) = mu(lin_dofs[l]);
      // And $\lambda_l(p_l+3) = 0.5$, $\lambda_{(l+1) mod 3}(p_l+3) = 0.5$;
      // $l =1,2,3$. Hence we copy 0.5 of the respective coefficients!
      zeta(quad_dofs[l + 3]) =
          0.5 * (mu(lin_dofs[l]) + mu(lin_dofs[(l + 1) % 3]));
    }
    /* END_SOLUTION */
  }
  return zeta;
}
/* SAM_LISTING_END_5 */

}  // namespace HandlingDOFs
