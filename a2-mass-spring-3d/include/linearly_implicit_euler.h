#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <EigenTypes.h>

template <typename FORCE, typename STIFFNESS>
inline void linearly_implicit_euler(
    Eigen::VectorXd &q,
    Eigen::VectorXd &qdot,
    double dt,
    const Eigen::SparseMatrixd &mass,
    FORCE &force,
    STIFFNESS &stiffness,
    Eigen::VectorXd &tmp_force,
    Eigen::SparseMatrixd &tmp_stiffness)
{
    // todo
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    stiffness(tmp_stiffness, q, qdot);
    force(tmp_force, q, qdot);
    solver.compute(mass - dt * dt * tmp_stiffness);
    qdot = solver.solve(mass * qdot + dt * tmp_force);
    q = q + dt * qdot;
}
