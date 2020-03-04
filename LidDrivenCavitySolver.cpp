#include <iostream>
using namespace std;

#include "LidDrivenCavity.h"
// #include "cxxopts.hpp"

int main(int argc, char **argv)
{
    // Create instance of options
    // cxxopts::Options options("Lid-Driven Cavity", "Vorticity StreamFunction");
    // options.add_options();

    // Create a new instance of the LidDrivenCavity class
    // Uses a pointer so solver is actually storing the address of the new instance
    
    // Check C++ version
    
    double Lx_arg = stod(argv[1]);
    double Ly_arg = stod(argv[2]);
    int Nx_arg = atoi(argv[3]);
    int Ny_arg = atoi(argv[4]);
    int Px_arg = atoi(argv[5]);
    int Py_arg = atoi(argv[6]);
    double dt_arg = stod(argv[7]);
    double T_arg = stod(argv[8]);
    double Re_arg = stod(argv[9]);

    LidDrivenCavity* solver = new LidDrivenCavity();
    
    // Configure the solver here...
    // Arrow dereferences and initialises the solver at that address
    solver -> Initialise(Lx_arg, Ly_arg, Nx_arg, Ny_arg, Px_arg, Py_arg, dt_arg, T_arg, Re_arg);

    // Run the solver
    solver -> Integrate();

	return 0;
}