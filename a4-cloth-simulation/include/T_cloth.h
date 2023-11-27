#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/**
 * Input:
 * @param: qdot - generalized velocities for the FEM system
 * @param: V - the nx3 matrix of undeformed vertex positions
 * @param: F - the mx3 matrix of triangle-vertex indices
 * @param: M - the sparse mass matrix for the FEM system.
 * Output:
 * @param: T - the kinetic energy for the entire cloth mesh.
 */
void T_cloth(
    double &T,
    Eigen::Ref<const Eigen::VectorXd> qdot,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::MatrixXi> F,
    Eigen::SparseMatrixd &M);