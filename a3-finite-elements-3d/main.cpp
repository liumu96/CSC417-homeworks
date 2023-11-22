#include <iostream>
#include <thread>

#include "assignment_setup.h"
#include "visualization.h"

// Simulate State
Eigen::VectorXd q;
Eigen::VectorXd qdot;

// simulation time and time step
double t = 0;     // simulation time
double dt = 0.01; // time step

// simulation loop
bool simulating = true;

bool simulation_callback();
bool draw_callback(igl::opengl::glfw::Viewer &viewer);

int main(int argc, char **argv)
{
    assignment_setup(argc, argv, q, qdot);

    // run simulation thread to avoid slowing down the Ui
    std::thread simulation_thread(simulation_callback);
    simulation_thread.detach();

    // setup libigl viewer and activate
    Visualize::setup(q, qdot, true);
    Visualize::viewer().callback_post_draw = &draw_callback;
    Visualize::viewer().launch();

    std::cout << "Finite Elements" << std::endl;
    return 0;
}

bool draw_callback(igl::opengl::glfw::Viewer &viewer)
{
    draw(q, qdot, t);
    return false;
}

bool simulation_callback()
{
    while (simulating)
    {
        simulate(q, qdot, dt, t);
        t += dt;
    }
    return false;
}
