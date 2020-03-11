// Header Files
#include "LidDrivenCavity.h"
#include "PoissonSolver.h"

// Basic
#include <iostream>
#include <exception>
#include <algorithm>
#include <math.h>
#include <fstream>

// External Libraries
#include "cblas.h"

using namespace std;

// Constructor (will be left empty)
LidDrivenCavity::LidDrivenCavity()
{
}

// Destructor (will be left empty)
LidDrivenCavity::~LidDrivenCavity()
{
}

// Setting the cavity size
void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    // Catch and return errors if either value is not positive
    if (xlen <= 0) {
        throw domain_error("Lx must be a positive value. Program terminated");
    }
    if (ylen <= 0) {
        throw domain_error("Ly must be a positive value. Program terminated");
    }

    // Set values if checks pass
    Lx = xlen;
    Ly = ylen;
}

// Setting discretised grid size
void LidDrivenCavity::SetGridSize(double nx, double ny)
{
    // Catch and return erros if either value is not a positive integer
    if (nx - int(nx) != 0 || nx <= 0) {
        throw domain_error("Nx must be a positive integer. Program terminated");
    }
    if (ny - int(ny) != 0 || ny <= 0) {
        throw domain_error("Ny must be a positive integer. Program terminated");
    }

    // Set values if checks pass
    Nx = int(nx);
    Ny = int(ny);
}

// Setting final time
void LidDrivenCavity::SetFinalTime(double finalt)
{
    // Catch and return error if time input is not positive
    if (finalt <= 0) {
        throw domain_error("T must be a positive value. Program terminated");
    }

    // Set value if check passes
    T = finalt;
}

// Setting Reynolds number
void LidDrivenCavity::SetReynoldsNumber(double re)
{
    // Catch and return error if Reynolds number is not positive
    if (re <= 0) {
        throw domain_error("Re must be a positive value. Program terminated");
    }

    // Set value if check passes
    Re = re;
}

// Setting time step size
void LidDrivenCavity::SetTimeStep(double deltat, double Lx, double Ly, unsigned int Nx, unsigned int Ny, double Re)
{
    // Catch and return error if timestep is not positive
    if (deltat <= 0 || deltat >= Re*(Lx/Nx)*(Ly/Ny)/4) {
        throw domain_error("dt must be a positive value and within the allowed range. Program terminated");
    }

    // Set value if check passes
    dt = deltat;
}

// Setting number of partitions
void LidDrivenCavity::SetPartitions(double px, double py, int nx, int ny)
{
    // Catch and return erros if either value is not a positive integer, or fails to properly divide domain
    if (px - int(px) != 0 || px <= 0 || (nx - 1)/px != int((nx - 1)/px)) {
        throw domain_error("Px must be a positive integer and appropriately divide the domain. Program terminated");
    }
    if (py - int(py) != 0 || py <= 0 || (ny - 1)/py != int((ny - 1)/py)) {
        throw domain_error("Py must be a positive integer and appropriately divide the domain. Program terminated");
    }

    // Set values if checks pass
    Px = int(px);
    Py = int(py);
}

// Initialising variables and outputting errors
void LidDrivenCavity::Solve(const double Lx_arg, const double Ly_arg, const double Nx_arg, const double Ny_arg,
                                const double Px_arg, const double Py_arg, const double dt_arg, const double T_arg, 
                                const double Re_arg)
{  
    // Perform argument checks and assign into variable if passed
    //--------------------------------------------------------------------------------------------------------------

    // Set cavity/domain size
    try {
        SetDomainSize(Lx_arg, Ly_arg);
    }
    catch (domain_error& e) {
        cout << "Error: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Set grid size
    try {
        SetGridSize(Nx_arg, Ny_arg);
    }
    catch (domain_error& e) {
        cout << "Error: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Set time step
    try {
        SetTimeStep(dt_arg, Lx_arg, Ly_arg, Nx_arg, Ny_arg, Re_arg);
    }
    catch (domain_error& e) {
        cout << "Error: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Set final time
    try {
        SetFinalTime(T_arg);
    }
    catch (domain_error& e) {
        cout << "Error: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Set Reynolds number
    try {
        SetReynoldsNumber(Re_arg);
    }
    catch (domain_error& e) {
        cout << "Error: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Set number of partitions
    try {
        SetPartitions(Px_arg, Py_arg, Nx_arg, Ny_arg);
    }
    catch (domain_error& e) {
        cout << "Error: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    //--------------------------------------------------------------------------------------------------------------

    // Implementing solver
    //--------------------------------------------------------------------------------------------------------------

    // Initial conditions
    const int Nx_const = Nx;
    const int Ny_const = Ny;

    double psi[Nx_const][Ny_const];
    double omega[Nx_const][Ny_const];
    double omega_new[Nx_const][Ny_const];

    // Psi and Omega initial zero matrix
    for (int i=0; i<Nx; i++){
        for (int j=0; j<Ny; j++) {
            psi[i][j] = 0;
            omega[i][j] = 0;
            omega_new[i][j] = 0;
        }
    }

    // Useful variables
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    double U = 1.0;
    int t_steps = ceil(T/dt);

    // Looping through every time increment
    for (int i=1; i<2; i++) {  // Change the max to t_steps when ready
        
        // Calculating vorticity boundary conditions at time t
        //---------------------------------------------------------------------------------------------------------
        for (int j=0; j<Nx; j++) {
            omega[0][j] = 2/(dy*dy) * (psi[0][j] - psi[1][j]) - 2*U/dy;  // Top Surface
            omega_new[0][j] = omega[0][j];
            omega[Nx-1][j] = 2/(dy*dy) * (psi[Nx-1][j] - psi[Nx-2][j]);  // Bottom Surface
            omega_new[Nx-1][j] = omega[Nx-1][j];
        }
        for (int j=1; j<Ny-1; j++) {
            omega[j][0] = 2/(dx*dx) * (psi[j][0] - psi[j][1]);           // Left Surface
            omega_new[j][0] = omega[j][0];
            omega[j][Nx-1] = 2/(dx*dx) * (psi[j][Nx-1] - psi[j][Nx-1]);  // Right Surface
            omega_new[j][Nx-1] = omega[j][Nx-1];
        }
        //---------------------------------------------------------------------------------------------------------

        // Calculate interior vorticity at time t
        //---------------------------------------------------------------------------------------------------------
        for (int j=1; j<Nx-1; j++) {
            for (int k=1; k<Ny-1; k++) {
                omega[j][k] = - ( (psi[j-1][k] - 2*psi[j][k] + psi[j+1][k])/(dx*dx) + 
                                  (psi[j][k+1] - 2*psi[j][k] + psi[j][k-1])/(dy*dy) );
            }
        }
        //---------------------------------------------------------------------------------------------------------

        // Calculate interior vorticity at time t + dt
        //---------------------------------------------------------------------------------------------------------
        for (int j=1; j<Nx-1; j++) {
            for (int k=1; k<Ny-1; k++) {
                omega_new[j][k] = ( (1/Re)*( (omega[j-1][k]-2*omega[j][k]+omega[j+1][k])/(dx*dx) + 
                                    (omega[j][k+1]-2*omega[j][k]+omega[j][k-1]) ) + 
                                    (psi[j-1][k]-psi[j+1][k])/(2*dx)*(omega[j][k+1]-omega[j][k-1])/(2*dy) -
                                    (psi[j][k+1]-psi[j][k-1])/(2*dy)*(omega[j-1][k]-omega[j+1][k])/(2*dx) ) * dt +
                                    omega[j][k];
            }
        }
        //---------------------------------------------------------------------------------------------------------

        // Solve the Poisson problem to calculate stream-function at time t + dt
        //---------------------------------------------------------------------------------------------------------        
        PoissonSolver* poisson = new PoissonSolver();
        poisson -> SolvePoisson((double*)omega_new, Ny, Nx);
        //---------------------------------------------------------------------------------------------------------
    }

    // Writing the vorticity and streamfunction values to a file for observation
    ofstream myfile1;
    ofstream myfile2;
    ofstream myfile3;
    myfile1.open("stream.txt");
    myfile2.open("vorticity.txt");
    myfile3.open("vorticity_new.txt");
    for (int i=0; i<Nx; i++){
        for (int j=0; j<Nx; j++) {
            myfile1 << psi[i][j] << " ";
            myfile2 << omega[i][j] << " ";
            myfile3 << omega_new[i][j] << " ";
        }
        myfile1 << endl;
        myfile2 << endl;
        myfile3 << endl;
    }
    myfile1.close();
    myfile2.close();
    myfile3.close();

}
