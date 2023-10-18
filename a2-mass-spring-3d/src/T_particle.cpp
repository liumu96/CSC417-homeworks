#include <T_particle.h>

// Compute the kinetic energy of a single mass particle
void T_particle(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, double mass)
{

    T = 1 / 2. * mass * qdot.dot(qdot);
}