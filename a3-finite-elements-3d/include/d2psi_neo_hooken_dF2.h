#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * @param psi - the 9x9 Hessian of the potential energy wrt to the deformation gradient
 * @param F - the dense 3x3 deformation gradient
 * @param C, D - material parameters for the Neo-Hookean model
 */
void d2psi_nei_hookean_dF2(
    Eigen::Matrix99d &ddw,
    Eigen::Ref<const Eigen::Matrix3d> F,
    double C,
    double D);