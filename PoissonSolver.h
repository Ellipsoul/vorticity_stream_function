#pragma once

#include <string>
using namespace std;

class PoissonSolver
{
    public:
    PoissonSolver();
    ~PoissonSolver();

    void SolvePoisson(double * omega_new, int Ny, int Nx, double dx, double dy);
    void PassPoisson(int Nx, int Ny, double* psi_new);

    private:
    double* b;
    int n;
};