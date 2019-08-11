/* **********************************************************************
   Demo code for course "Numerical Methods for PDEs
   Section "Case Study: Triangular Linear FEM in Two Dimensions
   ********************************************************************** */
#ifndef SLFEM2D
#define SLFEM2D

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>

#include <exception>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

/**
 * @brief The TriaMesh2D struct describes the triangulation in the form of a
 * Coordinates matrix that holds the coordinates of vertex i in its i-th row and
 * an Elements matrix that holds the three indices of the vertices forming a
 * triangle in each row
 */
/* SAM_LISTING_BEGIN_1 */
// Matrix containing vertex coordinates of a triangle
using TriGeo_t = Eigen::Matrix<double, 2, 3>;
struct TriaMesh2D {
  // Constructor: reads mesh data from file, whose name is passed
  TriaMesh2D(const std::string&);  // \Label[line]{tm:cs}
  virtual ~TriaMesh2D(void) {}

  // Retrieve coordinates of vertices of a triangles as rows of a 3x2 matrix
  TriGeo_t getVtCoords(std::size_t) const;

  // Data members describing geometry and topolgy
  Eigen::Matrix<double, Eigen::Dynamic, 2> Coordinates;
  Eigen::Matrix<int, Eigen::Dynamic, 3> Elements;
};
/* SAM_LISTING_END_1 */

/**
 * @brief The MatrixAssembler class assembles a Galerkin matrix given a
 * function for computing local matrices from TriGeo_t triangles.
 */

typedef std::function<Eigen::Matrix3d(const TriGeo_t&)> LocalMatrixHandle_t;

typedef std::vector<Eigen::Triplet<double>> Triplet_t;

class MatrixAssembler {
 public:
  // Constructor: stores element matrix assembler
  MatrixAssembler(LocalMatrixHandle_t getElementMatrix)
      : localMatrixHandle(std::move(getElementMatrix)){};

  // Assemble the Galerkin matrix for the provided mesh
  Eigen::SparseMatrix<double> Assemble(TriaMesh2D const& mesh);

 private:
  LocalMatrixHandle_t localMatrixHandle;
};

typedef std::function<Eigen::Matrix3d(const TriGeo_t&)> LocalMatrixHandle_t;

typedef std::vector<Eigen::Triplet<double>> Triplet_t;

/**
 * @brief The SlowMatrixAssembler class assembles a Galerkin matrix given a
 * function for computing local matrices from TriGeo_t triangles using a less
 * efficient way to build the sparse matrix.
 */
class SlowMatrixAssembler {
 public:
  // Constructor: stores element matrix assembler
  SlowMatrixAssembler(LocalMatrixHandle_t getElementMatrix)
      : localMatrixHandle(std::move(getElementMatrix)){};

  // Assemble the Galerkin matrix for the provided mesh
  Eigen::SparseMatrix<double> Assemble(TriaMesh2D const& mesh);

 private:
  LocalMatrixHandle_t localMatrixHandle;
};

typedef std::function<double(const Eigen::Vector2d&)> FHandle_t;

typedef std::function<Eigen::Vector3d(const TriGeo_t&, FHandle_t)>
    LocalVectorHandle_t;

/**
 * @brief The VectorAssembler class assembles a load vector given a
 * function for computing local vectors from TriGeo_t triangles.
 */
class VectorAssembler {
 public:
  // Constructor: stores element matrix assembler and load function
  VectorAssembler(LocalVectorHandle_t const& getElementVector,
                  const FHandle_t& sourceFunction)
      : localVectorHandle(std::move(getElementVector)),
        sourceFunction(sourceFunction){};

  // Assemble the load vector for the provided mesh
  Eigen::VectorXd Assemble(TriaMesh2D const& mesh);

 private:
  LocalVectorHandle_t localVectorHandle;
  const FHandle_t& sourceFunction;
};

/**
 * @brief The finite element solver solves the problem
 * \int_{\Omega} grad(u)grad(v) + uv dx = \int_{\Omega} fv dx
 * using MatrixAssembler and VectorAssembler. If the inefficient_flag is
 * set SlowMatrixAssembler is used to demonstrate the negative effect of not
 * building the Galerkin matrix from triplets.
 */
class FESolver {
 public:
  // Constructor: stores source function f and a flag indicating
  // if the slow Galerkin matrix assembler should be used
  FESolver(const FHandle_t& sourceFunction, int inefficient_flag = 0)
      : sourceFunction(sourceFunction) {
    inefficient = inefficient_flag;
  };

  // Solve the discretized system
  Eigen::VectorXd Solve(TriaMesh2D const& mesh);

 private:
  const FHandle_t& sourceFunction;
  int inefficient;
};
#endif
