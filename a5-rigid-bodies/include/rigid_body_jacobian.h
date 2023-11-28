#pragma once

#include <Eigen/Dense>
#include "EigenTypes.h"

/**
 * Input:
 * @param: R - rotation matrix for rigid body
 * @param: p - world space position of center-of-mass
 * @param: X -  undeformed position at which to compute the Jacobian.
 * Output:
 * @param: J - the rigid body jacobian acting on the undeformed space point X.
 */
void rigid_body_jacobian(
    Eigen::Matrix36d &J,
    Eigen::Ref<const Eigen::Matrix3d> R,
    Eigen::Ref<const Eigen::Vector3d> p,
    Eigen::Ref<const Eigen::Vector3d> X);