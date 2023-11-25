#include "dV_spring_particle_particle_dq.h"

void dV_spring_particle_particle_dq(
    Eigen::Ref<Eigen::Vector6d> f,
    Eigen::Ref<const Eigen::Vector3d> q0,
    Eigen::Ref<const Eigen::Vector3d> q1,
    double l0,
    double stiffness)
{
    f.setZero();
    // compute the spring direction（q1 - q0）
    Eigen::Vector3d delta_q = q1 - q0;

    // compute the spring length ||q1 - q0||
    double norm_delta_q = delta_q.norm();

    // compute the magnitude of force k * (||q1 - 10|| - l0)
    double force_magnitude = stiffness * (norm_delta_q - l0);

    // Calculate the component of the force on each node
    Eigen::Vector3d force = force_magnitude * (delta_q / norm_delta_q);

    // construct the force
    f.head(3) = -force; // 对 q0 的梯度
    f.tail(3) = force;  // 对 q1 的梯度
}