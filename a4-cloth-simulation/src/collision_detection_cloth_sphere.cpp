#include "collision_detection_cloth_sphere.h"

// Detect if any mesh vertex falls inside a sphere centered at (0,0,0.4) with radius 0.22
void collision_detection_cloth_sphere(
    std::vector<unsigned int> &cloth_index,
    std::vector<Eigen::Vector3d> &normals,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::Vector3d> center,
    double radius)
{
    cloth_index.clear();
    normals.clear();

    for (int i = 0; i < q.size() / 3; i++)
    {
        Eigen::Vector3d x = q.segment<3>(i * 3);
        if ((x - center).norm() <= radius)
        {
            // Collision happens!
            cloth_index.push_back(i);
            normals.push_back((x - center).normalized());
        }
    }
}