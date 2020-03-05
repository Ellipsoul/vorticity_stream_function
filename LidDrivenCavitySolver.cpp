// Basic
#include <iostream>

// External Libraries and Files
#include "LidDrivenCavity.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
    // Create a new instance of the LidDrivenCavity class
    // Uses a pointer so solver is actually storing the address of the new instance
    LidDrivenCavity* solver = new LidDrivenCavity();
    
    // Specify options that will be available to user

    // Option description
    po::options_description opts(
        "Numerically solves the lid-driven cavity problem using the vorticity stream function formulation"
    );
    // Options (no default values)
    opts.add_options()
        ("Lx", po::value<double>(),       "Length of the domain in the x-direction")
        ("Ly", po::value<double>(),       "Length of the domain in the y-direction")
        ("Nx", po::value<unsigned int>(), "Number of grid points in x-direction")
        ("Ny", po::value<unsigned int>(), "Number of grid points in y-direction")
        ("Px", po::value<unsigned int>(), "Number of partitions in the x-direction (parallel)")
        ("Py", po::value<unsigned int>(), "Number of partitions in the y-direction (parallel)")
        ("dt", po::value<double>(),        "Time step size")
        ("T",  po::value<double>(),        "Final time")
        ("Re", po::value<double>(),        "Reynolds number");

    // Parse command line options and generate map of options and values
    //po::variables_map vm;
    //po::store(po::parse_command_line(argc, argv, opts), vm);
    //po::notify(vm);

    // Extract parameters into the proper data type
    //const double Lx = vm["Lx"].as<double>();
    //const double Ly = vm["Ly"].as<double>();

    //std::cout << Lx << endl;


    // Configure the solver here...
    // Arrow dereferences and initialises the solver at that address
    solver -> Initialise();

    // Run the solver
    solver -> Integrate();

	return 0;
}