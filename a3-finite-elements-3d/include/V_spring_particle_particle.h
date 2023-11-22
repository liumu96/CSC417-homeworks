#pragma once

#include <Eigen/Dense>

/**
 * @param q0 - the generalized coordinates of the first node of the spring
 * @param q1 - the generalized coordinates of the second node of the spring
 * @param l0 - the undeformed length of the spring
 * @param stiffness - the stiffness constant for this spring
 * @param V - the potential energy of this spring
 */
void V_spring_particle_particle(
    double &V,
    Eigen::Ref<const Eigen::Vector3d> q0,
    Eigen::Ref<const Eigen::Vector3d> q1,
    double l0,
    double stiffness);