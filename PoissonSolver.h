#pragma once

#include <string>
using namespace std;

class PoissonSolver
{
    public:
    PoissonSolver();
    ~PoissonSolver();

    void SolvePoisson(double * omega_new, int Ny, int Nx);
};