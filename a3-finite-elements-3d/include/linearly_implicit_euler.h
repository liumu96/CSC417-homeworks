#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <EigenTypes.h>

/**
 * Input:
 * @param q - generalized coordinates for the FEM system
 * @param qdot - generalized velocity for the FEM system
 * @param dt - the time step in seconds
 * @param mass - the mass matrix
 * @param force(f, q, qdot) - - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
 * @param stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.
 * @param tmp_force - scratch space to collect forces
 * @param tmp_stiffness - scratch space to collect stiffness matrix
 * Output:
 * @param q - set q to the updated generalized coordinate using linearly implicit time integration
 * @param qdot - set q to the updated generalized coordinate using linearly implicit time integration
 */
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
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    stiffness(tmp_stiffness, q, qdot);
    force(tmp_force, q, qdot);
    solver.compute(mass - dt * dt * tmp_stiffness);
    qdot = solver.solve(mass * qdot + dt * tmp_force);
    q = q + dt * qdot;
}