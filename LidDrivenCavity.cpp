// Header File
#include "LidDrivenCavity.h"

// Basic
#include <iostream>
#include <exception>
#include <algorithm>

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

// Setting time step size
void LidDrivenCavity::SetTimeStep(double deltat)
{
    // Catch and return error if timestep is not positive
    if (deltat <= 0) {
        throw domain_error("dt must be a positive value. Program terminated");
    }

    // Set value if check passes
    dt = deltat;
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

void LidDrivenCavity::Initialise(const double Lx_arg, const double Ly_arg, const double Nx_arg, const double Ny_arg,
                                const double Px_arg, const double Py_arg, const double dt_arg, const double T_arg, 
                                const double Re_arg)
{  

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
        SetTimeStep(dt_arg);
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
}

void LidDrivenCavity::Integrate()
{
    // Initial Conditions
    
    // Initiliase streamfunction matrix of zeroes

}

