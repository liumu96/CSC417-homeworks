#pragma once

#include <Eigen/Dense>

template <typename FORCE>
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass, FORCE &force)
{
    Eigen::VectorXd f1, f2, f3, f4;
    Eigen::VectorXd q1 = q;
    Eigen::VectorXd qdot1 = qdot;
    force(f1, q1, qdot1);
    Eigen::VectorXd q2 = q + dt * qdot1 / 2;
    Eigen::VectorXd qdot2 = qdot + dt * f1 / mass / 2;
    force(f2, q2, qdot2);
    Eigen::VectorXd q3 = q + dt * qdot2 / 2;
    Eigen::VectorXd qdot3 = qdot + dt * f2 / mass / 2;
    force(f3, q3, qdot3);
    Eigen::VectorXd q4 = q + dt * qdot3;
    Eigen::VectorXd qdot4 = qdot + dt * f3 / mass;
    force(f4, q4, qdot4);
    q = q + dt * (qdot1 + 2 * qdot2 + 2 * qdot3 + qdot4) / 6;
    qdot = dt * (f1 + 2 * f2 + 2 * f3 + f4) / 6 / mass + qdot;
}