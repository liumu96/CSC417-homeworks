#pragma once

#include <Eigen/Dense>

template <typename FORCE>
inline void symplectic_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass, FORCE &force)
{
    // compute the force f
    Eigen::VectorXd f;
    force(f, q, qdot);

    qdot = qdot + dt * f / mass;
    q = q + dt * qdot;
}