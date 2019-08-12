/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "MyLinearLoadVector.h"

namespace ElementMatrixComputation {

MyLinearLoadVector::ElemVec MyLinearLoadVector::Eval(
    const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  const lf::base::size_type num_nodes{ref_el.NumNodes()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();

  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());

  return computeLoadVector(vertices, f_);
}

MyLinearLoadVector::ElemVec computeLoadVector(
    Eigen::MatrixXd vertices, const MyLinearLoadVector::FHandle_t f) {
  // Number of nodes of the element: triangles = 3, rectangles = 4
  const int num_nodes = vertices.cols();

  // Vector for returning element vector
  MyLinearLoadVector::elem_vec_t elem_vec = MyLinearLoadVector::elem_vec_t::Zero();

  /* BEGIN_SOLUTION */
  Eigen::MatrixXd midpoints(2, num_nodes);  // Matrix which will hold midpoint coordinates for triangle or quad.
  double area;  // Needed for quad rule.
  switch(num_nodes) {
    case 3 : {  // triangle
      area = 0.5 * ((vertices(0, 1) - vertices(0, 0)) * (vertices(1, 2) - vertices(1, 0)) -
                    (vertices(1, 1) - vertices(1, 0)) * (vertices(0, 2) - vertices(0, 0)));
      // midpoints is gonna ba 2 x 3 matrix
      auto a1 = vertices.col(0);
      auto a2 = vertices.col(1);
      auto a3 = vertices.col(2);
      midpoints << (a1(0) + a2(0)) / 2, (a2(0) + a3(0)) / 2, (a3(0) + a1(0)) / 2,
                   (a1(1) + a2(1)) / 2, (a2(1) + a3(1)) / 2, (a3(1) + a1(1)) / 2;
      break;;
    }
    case 4 : {  // quad
      area = (vertices(0, 1) - vertices(0, 0)) * (vertices(1, 3) - vertices(1, 0));
      // midpoints is gonna ba 2 x 4 matrix
      auto a1 = vertices.col(0);
      auto a2 = vertices.col(1);
      auto a3 = vertices.col(2);
      auto a4 = vertices.col(3);
      midpoints << (a1(0) + a2(0)) / 2, (a2(0) + a3(0)) / 2, (a3(0) + a4(0)) / 2, (a4(0) + a1(0)) / 2,
                   (a1(1) + a2(1)) / 2, (a2(1) + a3(1)) / 2, (a3(1) + a1(1)) / 2, (a4(1) + a1(1)) / 2;
      break;
    }
    default: {
      LF_ASSERT_MSG(false, "Wrong number of nodes for element.");
      break;
    }
  }

  // Evaluate functor f at midpoints
  Eigen::VectorXd f_midpoint_vals = Eigen::VectorXd::Zero(4);
  for (int i = 0; i < num_nodes; i++) {
    auto point = midpoints.col(i);
    f_midpoint_vals(i) = f(point);
  }

  // Assembly of element vector
  for (int i = 0; i < num_nodes; i++) {
    elem_vec(i) += 0.5 * f_midpoint_vals(i);
    elem_vec(i) += 0.5 * f_midpoint_vals((i + 2) % num_nodes);
  }

  // Multiply with factor acc. to Equ. 2.8.10
  double pre_factor = area / num_nodes;
  elem_vec *= pre_factor;


  /* END_SOLUTION */

  return elem_vec;
}

}  // namespace ElementMatrixComputation
