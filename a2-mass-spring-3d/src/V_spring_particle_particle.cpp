#include "V_spring_particle_particle.h"

void V_spring_particle_particle(
    double &V,
    Eigen::Ref<const Eigen::Vector3d> q0,
    Eigen::Ref<const Eigen::Vector3d> q1,
    double l0,
    double stiffness)
{
    // V = 1/2 * k(||q1 - 10|| - l0)^2
    double norm_diff = (q1 - q0).norm() - l0;
    V = 0.5 * stiffness * norm_diff * norm_diff;
}