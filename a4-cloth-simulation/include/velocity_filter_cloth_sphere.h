#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * Input:
 * @param: qdot - the 3nx1 generalized velocities of the cloth mesh
 * @param: index - a list of collision vertex indices from the collision detector
 * @param: normals - a list of collision normals from the collision detector
 * Output:
 * @param: qdot- the filtered 3nx1 generalized velocities
 */
void velocity_filter_cloth_sphere(
    Eigen::VectorXd &qdot,
    const std::vector<unsigned int> &index,
    const std::vector<Eigen::Vector3d> &normals);