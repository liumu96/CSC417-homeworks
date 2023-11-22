#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * Input:
 * @param V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
 * @param element - the 1x4 vertex indices for this tetrahedron
 * @param X - the position in the underformed space at which to compute the energy density
 * Output:
 * @param phi - the 4x1 values the basis functions
 */
void phi_linear_tetrahedron(
    Eigen::Vector4d &phi,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::RowVectorXi> element,
    Eigen::Ref<const Eigen::Vector3d> x);