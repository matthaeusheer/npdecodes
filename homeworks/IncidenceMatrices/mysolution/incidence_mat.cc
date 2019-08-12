#include "incidence_mat.h"

using namespace Eigen;
using lf::mesh::Mesh;

namespace IncidenceMatrices {

// @brief Create the mesh consisting of a triangle and quadrilateral
//        from the exercise sheet.
// @return Shared pointer to the hybrid2d mesh.
std::shared_ptr<Mesh> createDemoMesh() {
  // data type for a hybrid mesh in a world of dimension 2
  lf::mesh::hybrid2d::MeshFactory builder(2);

  // Add points
  builder.AddPoint(Vector2d{0, 0});    // (0)
  builder.AddPoint(Vector2d{1, 0});    // (1)
  builder.AddPoint(Vector2d{1, 1});    // (2)
  builder.AddPoint(Vector2d{0, 1});    // (3)
  builder.AddPoint(Vector2d{0.5, 1});  // (4)

  // Add the triangle
  // First set the coordinates of its nodes:
  MatrixXd nodesOfTria(2, 3);
  nodesOfTria << 1, 1, 0.5, 0, 1, 1;
  builder.AddEntity(
      lf::base::RefEl::kTria(),  // we want a triangle
      {1, 2, 4},                 // indices of the nodes
      std::make_unique<lf::geometry::TriaO1>(nodesOfTria));  // node coords

  // Add the quadrilateral
  MatrixXd nodesOfQuad(2, 4);
  nodesOfQuad << 0, 1, 0.5, 0, 0, 0, 1, 1;
  builder.AddEntity(lf::base::RefEl::kQuad(), {0, 1, 4, 3},
                    std::make_unique<lf::geometry::QuadO1>(nodesOfQuad));

  std::shared_ptr<Mesh> demoMesh = builder.Build();

  return demoMesh;
}

// @brief Compute the edge-vertex incidence matrix G for a given mesh
// @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
//             such as lf::mesh::hybrid2d::Mesh)
// @return The edge-vertex incidence matrix as Eigen::SparseMatrix<int>
SparseMatrix<int> computeEdgeVertexIncidenceMatrix(const Mesh& mesh) {
  // Store edge-vertex incidence matrix here
  assert(mesh.DimMesh() == 2 && "Mesh not 2D!");
  SparseMatrix<int, RowMajor> G;

  /* BEGIN_SOLUTION */
  // In every row (corresponding to an edge (co-dim 1 entity)) reserve space for two (incident nodes) non-zero entries
  const size_t n_edges = mesh.NumEntities(1);
  const size_t n_nodes = mesh.NumEntities(2);

  G = SparseMatrix<int, RowMajor>(n_edges, n_nodes);
  G.reserve(Eigen::VectorXi::Constant(n_edges, 2));
  // Loop over all edges and check their endpoint nodes.
  for (const auto& edge : mesh.Entities(1)) {
    auto edge_idx = mesh.Index(edge);
    auto nodes = edge.SubEntities(1);
    auto start_node_idx = mesh.Index(nodes[0]);
    auto end_node_idx = mesh.Index(nodes[1]);
    G.coeffRef(edge_idx, start_node_idx) += 1;
    G.coeffRef(edge_idx, end_node_idx) += -1;
  }
  /* END_SOLUTION */

  return G;
}

// @brief Compute the cell-edge incidence matrix D for a given mesh
// @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
//             such as lf::mesh::hybrid2d::Mesh)
// @return The cell-edge incidence matrix as Eigen::SparseMatrix<int>
SparseMatrix<int> computeCellEdgeIncidenceMatrix(const Mesh& mesh) {
  // Store cell-edge incidence matrix here
  SparseMatrix<int, RowMajor> D;

  /* BEGIN_SOLUTION */
  size_t n_cells = mesh.NumEntities(0);
  size_t n_edges = mesh.NumEntities(1);

  // Initialize sparse matrix and reserve space for max 4 non-zero entries per row (3 for triangle or 4 for squad)
  D = SparseMatrix<int, RowMajor>(n_cells, n_edges);
  D.reserve(Eigen::VectorXi::Constant(n_cells, 4));  // overestimating reserve size is ok

  // Loop over all entities of co-dim 0, which are cells (triangles and squads)
  for (const auto& cell : mesh.Entities(0)) {
    auto cell_idx = mesh.Index(cell);

    // Relative orientation of cell with corresponding edges gives 1 for same orientation, -1 for opposite orientation
    auto rel_orientations = cell.RelativeOrientations();

    // Loop over adjacent edges and add orientation sign to cell-edge incidence matrix.
    auto edges = cell.SubEntities(1);

    int iter_idx = 0;
    for (const auto& edge : edges) {
      size_t edge_idx = mesh.Index(edge);
      D.coeffRef(cell_idx, edge_idx) += to_sign(rel_orientations[iter_idx]);
      iter_idx++;
    }
  }
  /* END_SOLUTION */

  return D;
}

// @brief For a given mesh test if the product of cell-edge and edge-vertex
//        incidence matrix is zero: D*G == 0?
// @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
//             such as lf::mesh::hybrid2d::Mesh)
// @return true, if the product is zero and false otherwise
bool testZeroIncidenceMatrixProduct(const Mesh& mesh) {
  bool isZero = false;

  /* BEGIN_SOLUTION */
  SparseMatrix<int> G = computeEdgeVertexIncidenceMatrix(mesh);
  SparseMatrix<int> D = computeCellEdgeIncidenceMatrix(mesh);
  auto result_matrix = D * G;
  int norm = result_matrix.norm();
  isZero = norm == 0;
  /* END_SOLUTION */

  return isZero;
}

}  // namespace IncidenceMatrices
