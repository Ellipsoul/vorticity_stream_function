#pragma once

#include <string>
using namespace std;

/**
 * @class LidDrivenCavity
 * @brief Solves the Lid-Driven Cavity problem. Calls PoissonSolver class to solve the Poisson Problem
 */ 
class LidDrivenCavity
{
public:
    /// Class constructor and destructor
    LidDrivenCavity();
    ~LidDrivenCavity();

    /// Verifies and sets the cavity domain
    void SetDomainSize(double xlen, double ylen);
    /// Verifies and sets the discretisation grid size
    void SetGridSize(double nx, double ny);
    /// Verifies and sets the simulation duration
    void SetFinalTime(double finalt);
    /// Verifies and sets the Reynolds Number
    void SetReynoldsNumber(double Re);
    /// Verifies and sets the time increment step
    void SetTimeStep(double deltat, double Lx, double Ly, unsigned int Nx, unsigned int Ny, double Re);
    /// Verifies and sets the processor partitions for parallel programming
    void SetPartitions(double px, double py, int nx, int ny);
    /// Verifies the number of processors is consistent with partitions
    void VerifyProcessors(double px, double py);

    /// Performs all the above class functions and prepares the solver
    void Verify(const double Lx_arg, const double Ly_arg, const double Nx_arg, const double Ny_arg,
                    const double Px_arg, const double Py_arg, const double dt_arg, const double T_arg, 
                    const double Re_arg);
    /// Executes the solver
    void Solve();

private:
    // MPI variables
    int mpirank;    /// Processor rank
    bool mpiroot;   /// Identifies root processor
    int np;         /// Number of processors

    // Problem parameters
    double dt;          /// Time increment step
    double T;           /// Simulation duration
    int    Nx;          /// Grid points in x-direction
    int    Ny;          /// Grid points in y-direction
    int    Px;          /// Number of partitions in x-direction
    int    Py;          /// Number of partitions in y-direction
    double Lx;          /// Horizontal dimension of cavity
    double Ly;          /// Vertical dimension of cavity
    double Re;          /// Reynolds number
    double* psi_new;    /// Streamfunction array
};
