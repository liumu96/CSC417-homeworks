#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "EigenTypes.h"
#include "dV_membrane_corotational_dq.h"

/**
 * Input:
 * @param: q - generalized coordinates for the FEM system
 * @param: qdot - generalized velocity for the FEM system
 * @param: dX - an mx9 matrix which stores the flattened dphi/dX matrices for each tetrahedron.
 *       Convert this values back to 3x3 matrices using the following code (NOTE YOU MUST USE THE TEMPORARY VARIABLE tmp_row):
 *       Eigen::Matrix<double, 1,9> tmp_row
 *       tmp_row = dX.row(ei); //ei is the triangle index.
 *       Eigen::Map<const Eigen::Matrix3d>(tmp_row.data())
 * @param: V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
 * @param: F - the mx3 triangle connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
 * @param: a0 - the mx1 vector of undeformed triangle areas
 * @param: mu,lambda - material parameters for the cloth material model
 * Output:
 * @param: f - the vector 3xn vector of forces acting on each node of the mass-spring system
 */
void assemble_forces(
    Eigen::VectorXd &f,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::MatrixXd> qdot,
    Eigen::Ref<const Eigen::MatrixXd> dX,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::MatrixXi> F,
    Eigen::Ref<const Eigen::VectorXd> a0,
    double mu,
    double lambda);
