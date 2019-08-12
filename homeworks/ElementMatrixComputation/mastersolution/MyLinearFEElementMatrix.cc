/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "MyLinearFEElementMatrix.h"

namespace ElementMatrixComputation {

/* SAM_LISTING_BEGIN_1 */
MyLinearFEElementMatrix::ElemMat MyLinearFEElementMatrix::Eval(
    const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());
  // Matrix for returning element matrix
  MyLinearFEElementMatrix::elem_mat_t elem_mat;
  /* BEGIN_SOLUTION */
  // Initialize matrix containing the mass part of the element matrix
  MyLinearFEElementMatrix::elem_mat_t mass_elem_mat;
  // Retrieve laplace part of element matrix from LehrFEM's built-in laplace
  // element matrix builder
  auto laplace_elem_mat = laplace_elmat_builder_.Eval(cell);
  // Computations differ depending on the type of the cell
  switch (ref_el) {
    case lf::base::RefEl::kTria(): {
      double area = 0.5 * ((vertices(0, 1) - vertices(0, 0)) * (vertices(1, 2) - vertices(1, 0)) -
                           (vertices(1, 1) - vertices(1, 0)) * (vertices(0, 2) - vertices(0, 0)));
      // clang-format off
      mass_elem_mat << 2.0, 1.0, 1.0, 0.0, 
	1.0, 2.0, 1.0, 0.0, 
	1.0, 1.0, 2.0, 0.0, 
	0.0, 0.0, 0.0, 0.0;
      // clang-format on
      mass_elem_mat *= area / 12.0;
      break;
    }
    case lf::base::RefEl::kQuad(): {
      double area =
          (vertices(0, 1) - vertices(0, 0)) * (vertices(1, 3) - vertices(1, 0));
      // clang-format off
      mass_elem_mat << 4.0, 2.0, 1.0, 2.0, 
	2.0, 4.0, 2.0, 1.0, 
	1.0, 2.0, 4.0, 2.0, 
	2.0, 1.0, 2.0, 4.0;
      // clang-format on
      mass_elem_mat *= area / 36.0;
      break;
    }
    default: { LF_ASSERT_MSG(false, "Illegal cell type"); }
  }  // end switch
  elem_mat = laplace_elem_mat + mass_elem_mat;
  /* END_SOLUTION */
  return elem_mat;
}
/* SAM_LISTING_END_1 */
}  // namespace ElementMatrixComputation
