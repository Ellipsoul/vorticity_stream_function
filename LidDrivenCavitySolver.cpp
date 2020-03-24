// Basic
#include <iostream>
#include <typeinfo>

// External Libraries and Files
#include "LidDrivenCavity.h"
#include <boost/program_options.hpp>
#include <cblas.h>
#include <mpi.h>

// Setting up BLACS
extern "C" {
    // CBlacs Declarations
    void Cblacs_pinfo(int*, int*);
    void Cblacs_get(int, int, int*);
    void Cblacs_gridinit(int*, const char*, int, int);
    void Cblacs_pcoord(int, int, int*, int*);
    void Cblacs_gridexit(int);
    void Cblacs_barrier(int, const char*);
    void Cdgerv2d(int, int, int, double*, int, int, int);
    void Cdgesd2d(int, int, int, double*, int, int, int);
 
    int numroc_(int*, int*, int*, int*, int*);
}

namespace po = boost::program_options;

int main(int argc, char* argv[])
{

    // Initialise MPI
    //---------------------------------------------------------------------------------------------------------

    // Initialise MPI environment
    // MPI_Init(&argc, &argv);

    // Get number of processes
    // int world_size;
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get rank of process
    // int mpirank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

    // True (1) if process is the root process 
    // bool mpiroot = (mpirank == 0);
 
    // Helping variables
    // int iZERO = 0;
    //---------------------------------------------------------------------------------------------------------

    // Create a new instance of the LidDrivenCavity class
    // Uses a pointer so solver is actually storing the address of the new instance
    LidDrivenCavity* solver = new LidDrivenCavity();

    // Specify options that will be available to user
    // Option description
    po::options_description opts(
        "Numerically solves the lid-driven cavity problem using the vorticity stream function formulation"
    );
    // Adding options (without default values)
    opts.add_options()
        ("Lx", po::value<double>(), "Length of the domain in the x-direction")
        ("Ly", po::value<double>(), "Length of the domain in the y-direction")

        // Positive integer checks required
        ("Nx", po::value<double>(), "Number of grid points in x-direction")
        ("Ny", po::value<double>(), "Number of grid points in y-direction")
        ("Px", po::value<double>(), "Number of partitions in the x-direction (parallel)")
        ("Py", po::value<double>(), "Number of partitions in the y-direction (parallel)")

        ("dt", po::value<double>(), "Time step size")
        ("T",  po::value<double>(), "Final time")
        ("Re", po::value<double>(), "Reynolds number");

    // Parse command line options and generate map of options and values
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    // Extract command line arguments into parameters
    const double Lx_arg = vm["Lx"].as<double>();
    const double Ly_arg = vm["Ly"].as<double>();
    const double Nx_arg = vm["Nx"].as<double>();
    const double Ny_arg = vm["Ny"].as<double>();
    const double Px_arg = vm["Px"].as<double>();
    const double Py_arg = vm["Py"].as<double>();
    const double dt_arg = vm["dt"].as<double>();
    const double T_arg  = vm["T"].as<double>();
    const double Re_arg = vm["Re"].as<double>();

    // Verify and set command line inputs
    solver -> Verify(Lx_arg, Ly_arg, Nx_arg, Ny_arg, Px_arg, Py_arg, dt_arg, T_arg, Re_arg);

    cout << "All command line input checks passed" << endl;

    solver -> Solve();
 
 
	// MPI_Finalize();
}