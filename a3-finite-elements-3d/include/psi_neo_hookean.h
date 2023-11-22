#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * Input:
 * @param F - the dense 3x3 deformation gradient
 * @param C,D - material parameters for the Neo-Hookean model
 * Output:
 * @param psi - the neohookean energy
 */
void psi_neo_hookean(
    double &psi,
    Eigen::Ref<const Eigen::Matrix3d> F,
    double C,
    double D);