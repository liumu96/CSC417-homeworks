#include "velocity_filter_cloth_sphere.h"

// Project out components of the per-vertex velocities which are in the positive direction of the contact normal
void velocity_filter_cloth_sphere(
    Eigen::VectorXd &qdot,
    const std::vector<unsigned int> &indices,
    const std::vector<Eigen::Vector3d> &normals)
{
    for (int i = 0; i < indices.size(); i++)
    {
        Eigen::Vector3d normal = normals[i];
        Eigen::Vector3d v = qdot.segment<3>(indices[i] * 3);
        if (normal.dot(v) < 0)
        {
            qdot.segment<3>(indices[i] * 3) -= normal.dot(v) * normal;
        }
    }
}