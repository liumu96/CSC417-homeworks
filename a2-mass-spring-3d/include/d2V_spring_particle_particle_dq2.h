#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * @param q0 - the generalized coordiantes of the first node of the spring
 * @param q1 - the generalized coordiantes of the second node of the spring
 * @param l0 - the undeformed length of the spring
 * @param stiffness - the stiffness constant for this spring
 * @param H - the 6x6 dense per spring hessian of the energy function
 */
void d2V_spring_particle_particle_dq2(
    Eigen::Ref<Eigen::Matrix66d> H,
    Eigen::Ref<const Eigen::Vector3d> q0,
    Eigen::Ref<const Eigen::Vector3d> q1,
    double l0,
    double stiffness);