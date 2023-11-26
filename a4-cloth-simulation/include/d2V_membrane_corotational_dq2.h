#pragma once

#include <Eigen/Dense>
#include "EigenTypes.h"
#include "dsvd.h"

/**
 * @param: q - generalized coordinates for the FEM system
 * @param: dX - the 3x3 matrix containing dphi/dX
 * @param: V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
 * @param: element - the vertex indices of this triangle
 * @param: area - the area of this triangle
 * @param: mu,lambda - material parameters for the cloth material model
 * Output:
 * @param: H - the per-triangle Hessian of the potential energy (the linear model described in the README).
 */
void d2V_membrane_corotational_dq2(
    Eigen::Matrix99d &H,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::Matrix3d> dX,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::RowVectorXi> element,
    double area,
    double mu,
    double lambda);