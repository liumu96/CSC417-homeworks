#pragma once

#include <Eigen/Dense>

template <typename FORCE>
inline void forward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass, FORCE &force)
{
    // compute f
    Eigen::VectorXd f;
    force(f, q, qdot);
    // x(n+1) = x(n) + dt * v(n)
    q = q + dt * qdot;
    // v(n+1) = v(n) + dt * a = v(n) + dt * f / mass
    qdot = qdot + dt * f / mass;
}