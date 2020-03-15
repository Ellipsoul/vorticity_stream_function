#pragma once

#include <string>
using namespace std;

class LidDrivenCavity
{
public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(double nx, double ny);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void SetTimeStep(double deltat, double Lx, double Ly, unsigned int Nx, unsigned int Ny, double Re);
    void SetPartitions(double px, double py, int nx, int ny);

    void Solve(const double Lx_arg, const double Ly_arg, const double Nx_arg, const double Ny_arg,
                    const double Px_arg, const double Py_arg, const double dt_arg, const double T_arg, 
                    const double Re_arg);

    // Add any other public functions

private:
    // All private variables can be accessed by any member in the function
    double* v = nullptr;
    double* s = nullptr;

    double dt;
    double T;
    int    Nx;
    int    Ny;
    int    Px;
    int    Py;
    double Lx;
    double Ly;
    double Re;

    double* psi_new[];

};
