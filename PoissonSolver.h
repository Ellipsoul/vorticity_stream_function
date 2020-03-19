#pragma once

#include <string>
using namespace std;

class PoissonSolver
{
    double* b;
    int n;

    public:
        PoissonSolver();
        ~PoissonSolver();

        void SolvePoisson(double * omega_new, int Ny, int Nx, double dx, double dy);
        double* ReturnStream();
};