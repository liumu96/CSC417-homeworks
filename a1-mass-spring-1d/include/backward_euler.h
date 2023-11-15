#pragma once

#include <Eigen/Dense>

template <typename FORCE, typename STIFFNESS>
inline void backward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass, FORCE &force, STIFFNESS &stiffness)
{
    // compute f
    Eigen::VectorXd f;
    Eigen::MatrixXd k;
    force(f, q, qdot);
    stiffness(k, q, qdot);

    qdot = (qdot + dt * f / mass) / (1 - dt * dt * k(0, 0) / mass);
    // x(n+1) = x(n) + dt * v(n+1)
    q = q + dt * qdot;
}