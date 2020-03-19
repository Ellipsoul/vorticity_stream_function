#pragma once

#include <string>
using namespace std;

class PoissonSolver
{
    public:
        PoissonSolver();
        ~PoissonSolver();

        void SolvePoisson(double * omega_new, int Ny, int Nx, double dx, double dy);
        void ReturnStream(double * psi_new, int Nx, int Ny);

    private:
        double* A;
        double* b;
        int n;
        int* piv;
};