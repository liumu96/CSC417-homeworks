#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/**
 * Input:
 * @param: q - generalized coordinates for the FEM system
 * @param: V - the nx3 matrix of undeformed vertex positions
 * @param: F - the mx3 matrix of triangle-vertex indices
 * @param: density - the density of the cloth material
 * @param: areas - the mx1 vector of undeformed triangle areas
 * Output:
 * @param: M - sparse mass matrix for the entire mesh
 */
void mass_matrix_mesh(
    Eigen::SparseMatrixd &M,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::MatrixXi> F,
    double density,
    Eigen::Ref<const Eigen::VectorXd> areas);