#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * @param: q - generalized coordinates for the FEM system
 * @param: dX - the 3x3 matrix containing dphi/dX
 * @param: V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
 * @param: element - the vertex indices of this triangle
 * @param: area - the area of this triangle
 * @param: mu,lambda - material parameters for the cloth material model
 * Output:
 * @param: energy- the per-triangle potential energy (the linear model described in the README).
 */
void V_membrane_corotational(
    double &energy,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::Matrix3d> dX,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::RowVectorXi> element,
    double area,
    double mu,
    double lambda);