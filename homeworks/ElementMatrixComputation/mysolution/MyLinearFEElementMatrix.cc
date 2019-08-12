/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "MyLinearFEElementMatrix.h"

namespace ElementMatrixComputation {

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

  // Mass element matrix for triangle and quad case
  MyLinearFEElementMatrix::elem_mat_t mass_elem_mat;
  double area;
  switch(geo_ptr->RefEl()) {
    case lf::base::RefEl::kTria() : {
      std::cout << "Triangle!" << std::endl;
      // Compute area
      area = 0.5 * ((vertices(0, 1) - vertices(0, 0)) * (vertices(1, 2) - vertices(1, 0)) -
                    (vertices(1, 1) - vertices(1, 0)) * (vertices(0, 2) - vertices(0, 0)));
      // Fill mass_elem_mat with const matrix for triangle
      mass_elem_mat << 2, 1, 1, 0,
                       1, 2, 1, 0,
                       1, 1, 2, 0,
                       0, 0, 0, 0;
      break;
    }

    case lf::base::RefEl::kQuad() : {
      std::cout << "Quad!" << std::endl;
      // Compute area
      area = (vertices(0, 1) - vertices(0, 0)) * (vertices(1, 3) - vertices(1, 0));
      // Fill mass_elem_mat with const matrix for quad
      mass_elem_mat << 4, 2, 1, 2,
                       2, 4, 2, 1,
                       1, 2, 4, 2,
                       2, 1, 2, 4;
      break;
    }

    default :
      LF_ASSERT_MSG(false, "Cell type not supported.")
  }
  // Multiply by area since this is what element matrix depends on
  mass_elem_mat *= area / 36;
  // The laplace part
  auto elem_mat_laplace = laplace_elmat_builder_.Eval(cell);
  // Add both contributions, from laplace term and from second, analytic term
  elem_mat = elem_mat_laplace + mass_elem_mat;
  /* END_SOLUTION */

  return elem_mat;
}
}  // namespace ElementMatrixComputation