#include "dV_spring_particle_particle_dq.h"

void dV_spring_particle_particle_dq(Eigen::VectorXd &dV, const Eigen::VectorXd &q, double stiffness)
{
    dV.resize(1);
    // V = 0.5 * k * q^2 -> dV/dq = k * q
    dV = stiffness * q;
}