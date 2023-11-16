#pragma once

#include <Eigen/Dense>

void dV_gravity_particle_dq(
    Eigen::Ref<Eigen::Vector3d> f,
    double mass,
    Eigen::Ref<const Eigen::Vector3d> g);