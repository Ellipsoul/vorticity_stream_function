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
        int mype;           /// Processor rank
        int npe;            /// Number of processors
        int ctx;            /// Context
        int nrow;           /// Number of processor rows
        int ncol;           /// Number of processor columns
        int myrow;          /// Current processor row
        int mycol;          /// Current processor column
        char order;         /// Defines order of processors

        int n;              /// Columns in banded A matrix (global)
        int kl;             /// Lower diagonal bandwidth (global)
        int ku;             /// Upper diagonal bandwidth (global)
        int ldab;           /// Rows in banded A matrix (global, SCALAPACK form)

        //double* vorticity_vec;

        double* A;          /// Declare A matrix (global)
        double* b;          /// Declare b array (global)

        double* A_loc;      /// Declare A matrix (local)
        double* x;          /// Declare x array (local)
        int* ipiv;          /// Pivoting array
        double* work;       /// Delcare workspace
};