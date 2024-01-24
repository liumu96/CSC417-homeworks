#include "dV_spring_particle_particle_dq.h"

void dV_spring_particle_particle_dq(
    Eigen::Ref<Eigen::Vector6d> f,
    Eigen::Ref<const Eigen::Vector3d> q0,
    Eigen::Ref<const Eigen::Vector3d> q1,
    double l0, double stiffness)
{
    f.setZero();
    double dV_constant = stiffness * (1 - l0 / (q1 - q0).norm());

    f << dV_constant * (q0 - q1),
        dV_constant * (q1 - q0);
}