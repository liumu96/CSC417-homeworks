
#include <iostream>

#include "visualization.h"

#include "dV_spring_particle_particle_dq.h"
#include "d2V_spring_particle_particle_dq2.h"
#include "forward_euler.h"
#include "backward_euler.h"
#include "symplectic_euler.h"
#include "runge_kutta.h"

#include <Eigen/Dense>
#include <igl/readOBJ.h>

// simulation parameters
Eigen::VectorXd q;
Eigen::VectorXd q_dot;

double mass = 1.0;
double stiffness = 100.0;
double dt = 1e-2;
int integrator_type = 0;

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods);
bool simulate(igl::opengl::glfw::Viewer &viewer);

int main(int argc, char *argv[])
{

    std::cout << "mass-spring-1d" << std::endl;

    if (argc > 1)
    {
        if (argv[1] == std::string("be"))
        {
            integrator_type = 1; // backward-euler
        }
        else if (argv[1] == std::string("se"))
        {
            integrator_type = 2; // simplectic-euler
        }
        else if (argv[1] == std::string("rk"))
        {
            integrator_type = 3; // runge-kutta
        }
    }

    // load data for animations
    Eigen::MatrixXd V_cow, V_spring;
    Eigen::MatrixXi F_cow, F_spring;

    igl::readOBJ("../data/spot.obj", V_cow, F_cow);
    igl::readOBJ("../data/spring.obj", V_spring, F_spring);

    // setup simulation variables
    q.resize(1);
    q_dot.resize(1);

    q_dot(0) = 0;
    q(0) = 1;

    // setup libigl viewer and activate
    Visualize::setup(q, q_dot);
    Visualize::viewer().callback_post_draw = &simulate;
    Visualize::viewer().callback_key_down = &key_down;
    Visualize::add_object_to_scene(V_cow, F_cow, Eigen::RowVector3d(244, 165, 130) / 255.);
    Visualize::add_object_to_scene(V_spring, F_spring, Eigen::RowVector3d(200, 200, 200) / 255.);
    Visualize::viewer().launch();

    return 0;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
    switch (key)
    {
    case ' ':
        viewer.core().is_animating = !viewer.core().is_animating;
        return true;
    }
    return false;
}

bool simulate(igl::opengl::glfw::Viewer &viewer)
{
    if (Visualize::viewer().core().is_animating)
    {
        // take a time step
        auto force = [](
                         Eigen::VectorXd &f,
                         Eigen::VectorXd &q,
                         Eigen::VectorXd &qdot)
        {
            dV_spring_particle_particle_dq(f, q, stiffness);
            f *= -1;
        };

        auto stiff = [](
                         Eigen::MatrixXd &k,
                         Eigen::VectorXd &q,
                         Eigen::VectorXd &qdot)
        {
            d2V_spring_particle_particle_dq2(k, q, stiffness);
            k *= -1;
        };

        switch (integrator_type)
        {
        case 0:
            forward_euler(q, q_dot, dt, mass, force);
            break;
        case 1:
            backward_euler(q, q_dot, dt, mass, force, stiff);
            break;
        case 2:
            symplectic_euler(q, q_dot, dt, mass, force);
            break;
        case 3:
            runge_kutta(q, q_dot, dt, mass, force);
            break;

        default:
            break;
        }
        // update mesh positions
        Visualize::rigid_transform_1d(0, q(0));
        Visualize::scale_x(1, q(0));
    }

    return false;
}