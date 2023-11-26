#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * @param: q - generalized coordinates for the FEM system
 * @param: center - the position of the sphere center in the world space
 * @param: radius - the radius of the sphere in the world space
 * Output:
 * @param: cloth_index - the indices of all vertices currently in contact with the sphere
 * @param: normals - the outward facing contact normals for each contacting vertex.
 */
void collision_detection_cloth_sphere(
    std::vector<unsigned int> &cloth_index,
    std::vector<Eigen::Vector3d> &normals,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::Vector3d> center,
    double radius);