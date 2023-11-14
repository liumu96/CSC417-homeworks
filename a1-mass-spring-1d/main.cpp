
#include <iostream>

// #include <Eigen/Dense>

// #include <igl/readOBJ.h>

int integrator_type = 0;

bool simulate()
{

    switch (integrator_type)
    {
    case 0:
        /* code */
        break;
    case 1:
        /* code */
        break;
    case 2:
        /* code */
        break;
    case 3:
        /* code */
        break;

    default:
        break;
    }
    return false;
}

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

    // setup simulation variables

    // setup libigl viewer and activate

    return 0;
}