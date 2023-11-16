#include "T_particle.h"

void T_particle(
    double &T,
    Eigen::Ref<const Eigen::VectorXd> qdot,
    double mass)
{
    // Compute the kinetic energy of a single mass particle: T = 1/2 * m * v^2
    T = mass * qdot.dot(qdot) / 2.;
}