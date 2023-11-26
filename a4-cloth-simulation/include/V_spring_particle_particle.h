#pragma once

#include <Eigen/Dense>

/**
 * Input:
 * @param: q0 - the generalized coordinates of the first node of the spring
 * @param: q1 - the generalized coordinates of the second node of the spring
 * @param: l0 - the undeformed length of the spring
 * @param: stiffness - the stiffness constant for this spring
 * Output:
 * @param: V - potential energy for the spring
 */
void V_spring_particle_particle(
    double &V,
    Eigen::Ref<const Eigen::Vector3d> q0,
    Eigen::Ref<const Eigen::Vector3d> q1,
    double l0,
    double stiffness);