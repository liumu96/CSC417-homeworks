#pragma once

#include <Eigen/Dense>

/**
 * @param q - generalized coordinate of a particle
 * @param mass - the mass of particles in the mass spring system
 * @param g - the gravity acceleration vector
 * @param V - the potential energy due to gravity acting on this particle
 */
void V_gravity_particle(
    double &V,
    Eigen::Ref<const Eigen::Vector3d> q,
    double mass,
    Eigen::Ref<const Eigen::Vector3d> g);