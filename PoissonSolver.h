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
        // MPI Variables
        int mpirank;
        bool mpiroot;

        // BLACS Variables
        int mype, npe, ctx, nrow, ncol, myrow, mycol;
        char order; 

        int n;         // Columns in banded A matrix (global)
        int kl;        // Low diagonal bandwidth (global)
        int ku;        // Upper diagonal bandwidth (global)
        int ldab;      // Rows in banded A matrix (global, SCALAPACK form)
        double* A;     // Declare A matrix (global)
        double* b;     // Declare b array (global)

        double* A_loc; // Declare A matrix (local)
        double* x;     // Declare x array (local)
        int* ipiv;     // Pivoting array
        double* work;  // Delcare workspace

};