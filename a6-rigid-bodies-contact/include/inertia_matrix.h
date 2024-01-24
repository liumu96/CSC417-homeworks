#pragma once

#include <Eigen/Dense>
#include "EigenTypes.h"

/**
 * Input:
 * @param: V - the nx3 matrix of vertices.
 * @param: F - the mx3 matrix of triangle vertex indices.
 * @param: density - the material density.
 * Output:
 * @param: I - the 3x3 angular inertia matrix
 * @param: center - the center of mass of the object
 * @param: mass - the total mass of the object
 */
void inertia_matrix(
    Eigen::Matrix3d &I,
    Eigen::Vector3d &center,
    double &mass,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::MatrixXi> F,
    double density);