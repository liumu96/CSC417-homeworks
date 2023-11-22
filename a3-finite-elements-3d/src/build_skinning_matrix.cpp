#include "build_skinning_matrix.h"
#include "phi_linear_tetrahedron.h"
#include <vector>
#include <iostream>
#include <igl/AABB.h>
#include <igl/in_element.h>

void build_skinning_matrix(
    Eigen::SparseMatrixd &N,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::MatrixXi> T,
    Eigen::Ref<const Eigen::MatrixXd> V_skin)
{
    // Reference https://libigl.github.io/tutorial/#point-location

    // Initialize an AABB tree for efficient point location within the tetrahedral mesh
    igl::AABB<Eigen::MatrixXd, 3> tree;
    tree.init((Eigen::MatrixXd)V, (Eigen::MatrixXi)T);
    // Vector to store the indices of the tetrahedra containing each vertex in V_skin
    Eigen::VectorXi I;

    // Use the AABB tree to efficiently find the containing tetrahedra for each vertex in V_skin
    igl::in_element((Eigen::MatrixXd)V, (Eigen::MatrixXi)T, (Eigen::MatrixXd)V_skin, tree, I);

    // Resize the sparse matrix N to store the skinning weights
    N.resize(V_skin.rows(), V.rows());
    N.setZero(); // Resize the sparse matrix N to store the skinning weights

    // Vector of triplets to store the non-zero entries of the sparse matrix
    std::vector<Eigen::Triplet<double>> N_entries;

    // Loop through each vertex in V_skin
    for (int i = 0; i < V_skin.rows(); i++)
    {
        // Vector to store the barycentric coordinates within the tetrahedron for the current V_skin vertex
        Eigen::Vector4d phi;

        // Get the indices of the tetrahedron containing the current V_skin vertex
        Eigen::RowVector4i bounding_tetrahedra = T.row(I(i));

        // Compute the barycentric coordinates using linear tetrahedral interpolation
        phi_linear_tetrahedron(phi, V, bounding_tetrahedra, V_skin.row(i).transpose());

        // Loop through the vertices of the bounding tetrahedron
        for (int vertex_i = 0; vertex_i < 4; vertex_i++)
        {
            // Add a triplet to the vector for each non-zero entry in the sparse matrix
            N_entries.push_back(Eigen::Triplet<double>(
                i,                             // Row index in the sparse matrix N
                bounding_tetrahedra(vertex_i), // Column index in the sparse matrix N
                phi(vertex_i)));               // Value at the specified row and column
        }
    }

    // Set the sparse matrix N using the collected triplets
    N.setFromTriplets(N_entries.begin(), N_entries.end());
}
