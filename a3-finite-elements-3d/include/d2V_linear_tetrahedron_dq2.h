#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * @param H - the 12x12 Hessian of the potential energy for a single tetrahedron
 * @param q - generalized coordinates for the FEM system
 * @param V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position
 * @param element - the 1x4 vertex indices for this tetrahedron
 * @param volume - the undeformed tetrahedron volume
 * @param C, D - material parameters for the Neo-Hookean model
 */
void d2V_linear_tetrahedron_dq2(
    Eigen::Matrix1212d &H,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::RowVectorXi> element,
    double volume,
    double C,
    double D);