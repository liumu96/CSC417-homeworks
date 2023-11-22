#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * @param M - dense 12x12 per-tetrahedron mass matrix
 * @param qdot - generalized velocity for the FEM system
 * @param T - the mx4 vertex indices for tet mesh
 * @param density - density of material
 * @param v0 - the undeformed tetrahedra volumes
 */
void mass_matrix_mesh(
    Eigen::SparseMatrixd &M,
    Eigen::Ref<const Eigen::VectorXd> qdot,
    Eigen::Ref<const Eigen::MatrixXi> T,
    double density,
    Eigen::Ref<const Eigen::VectorXd> v0);