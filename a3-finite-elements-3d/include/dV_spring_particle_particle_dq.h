#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * @param q0 - the generalized coordiantes of the first node of the spring
 * @param q1 - the generalized coordiantes of the second node of the spring
 * @param l0 - the undeformed length of the spring
 * @param stiffness - the stiffness constant for this spring
 * @param f - the 6x1 per spring generalized energy gradient vector
 */
void dV_spring_particle_particle_dq(
    Eigen::Ref<Eigen::Vector6d> f,
    Eigen::Ref<const Eigen::Vector3d> q0,
    Eigen::Ref<const Eigen::Vector3d> q1,
    double l0,
    double stiffness);