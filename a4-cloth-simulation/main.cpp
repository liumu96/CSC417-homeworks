#include <iostream>
#include <thread>

#include "assignment_setup.h"
#include "visualization.h"

// Simulate State
Eigen::VectorXd q;
Eigen::VectorXd qdot;

// Simulation time and time step
double t = 0;
double dt = 0.001;

// simulation loop
bool simulating = true;

bool simulation_callback();
bool draw_callback(igl::opengl::glfw::Viewer &viewer);

int main(int argc, char **argv)
{
    std::cout << "Start Cloth Simulation" << std::endl;

    // assignment specific setup
    assignment_setup(argc, argv, q, qdot);

    // run simulation in seperate thread to avoid slowing down the UI
    std::thread simulation_thread(simulation_callback);
    simulation_thread.detach();

    // setup libigl viewer and activate
    Visualize::setup(q, qdot, true);
    Visualize::viewer().callback_post_draw = *draw_callback;
    Visualize::viewer().launch();

    return 1;
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

bool draw_callback(igl::opengl::glfw::Viewer &viewer)
{
    draw(q, qdot, t);
    return false;
}
