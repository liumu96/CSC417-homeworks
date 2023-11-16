#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/**
 * @param qdot - generalized velocitu for the mass spring system
 * @param mass - the mass of a particle
 * @param T - kinetic energy of the all particles in the mass-spring system
 */
void T_particle(
    double &T,
    Eigen::Ref<const Eigen::VectorXd> qdot,
    double mass);