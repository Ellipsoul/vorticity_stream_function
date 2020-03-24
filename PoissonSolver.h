#pragma once

#include <string>
using namespace std;

/**
 * @class PoissonSolver
 * @brief Solves the Poisson problem and calculates the updated stream-function 
 */ 
class PoissonSolver
{
    public:
        /// Class constructor and destructor
        PoissonSolver();
        ~PoissonSolver();

        /// Executes the Poisson solver
        void SolvePoisson(double * omega_new, int Ny, int Nx, double dx, double dy);
        /// Returns the calculated stream-function
        void ReturnStream(double * psi_new, int Nx, int Ny);

    private:
        // MPI Variables
        int mpirank;        /// Processor rank
        bool mpiroot;       /// Identifies root processor

        // BLACS Variables
        int mype;
        int npe;
        int ctx;
        int nrow;
        int ncol;
        int myrow; 
        int mycol;
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