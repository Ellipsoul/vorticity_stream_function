#include "LidDrivenCavity.h"

#include <iostream>
using namespace std;

// Why does the class have a function too? What's it meant to do?
LidDrivenCavity::LidDrivenCavity()
{
}

// Destructor deleted the instance
LidDrivenCavity::~LidDrivenCavity()
{
}

// Difference between domain and grid?
void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
}

// Isn't deltat already the time step?
void LidDrivenCavity::SetTimeStep(double deltat)
{
}

// Same question
void LidDrivenCavity::SetFinalTime(double finalt)
{
}

// Same
void LidDrivenCavity::SetReynoldsNumber(double re)
{
}

void LidDrivenCavity::Initialise(double Lx_arg, double Ly_arg, int Nx_arg, int Ny_arg, int Px_arg, 
                                 int Py_arg, double dt_arg, double T_arg, double Re_arg)
{
    // These parameters can be accessed within this instance
    Lx = Lx_arg;
    Ly = Ly_arg;
    Nx = Nx_arg;
    Ny = Ny_arg;
    dt = dt_arg;
    T = T_arg;
    Re = Re_arg;

}

void LidDrivenCavity::Integrate()
{
    
}