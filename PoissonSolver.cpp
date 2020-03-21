// Header Files
#include "PoissonSolver.h"

// Basic
#include <iostream>
#include <math.h>
#include <fstream>
#include <cstdlib>

// External Libraries
#include <cblas.h>
#include <mpi.h>

using namespace std;

// Defining LAPACK routine for solving banded linear system
#define F77NAME(x) x##_
extern "C" {
    void F77NAME(pdgbsv) (const int& N, const int& BWL, const int& BWU, const int& NRHS, double * A,
                          const int& JA, int * desca, int * ipiv, double * x, const int& IB, int * descb,
                          double * work, const int& LW, int& info) ;
}

// Constructor
PoissonSolver::PoissonSolver()
{
}

// Destructor
PoissonSolver::~PoissonSolver()
{
    delete[] A;
    delete[] b;
    delete piv;
}

// Function to solve the Poisson problem
void PoissonSolver::SolvePoisson(double* omega_new, int Ny, int Nx, double dx, double dy) {
    
    // Visualise passed vorticity matrix into a file
    // if (MPI_ROOT) {
    //     ofstream myfile4;
    //     myfile4.open("vorticity_matrix_old_trans.txt");
    //     for (int i=0; i<Ny; i++) {
    //         for (int j=0; j<Nx; j++) {
    //             myfile4 << *(omega_new + i*Nx + j) << " ";
    //         }
    //         myfile4 << endl;
    //     }
    //     myfile4.close();
    // }

    // Parallel Banded Matrix Solver
    //---------------------------------------------------------------------------------------------------------
    int info;                                   // Status value
    const int N    = (Nx-2)*(Ny-2);             // Total problem size
    const int NB   = 4;                         // Blocking size (number of columns per process)
    const int BWL  = Nx-2;                      // Lower bandwidth
    const int BWU  = Nx-2;                      // Upper bandwidth
    const int NRHS = 1;                         // Number of RHS to solve
    const int JA   = 1;                         // Start offset in matrix (to use just a submatrix)
    const int IB   = 1;                         // Start offset in RHS vector (to use just a subvector)
    const int LA   = (1 + 2*BWL + 2*BWU)*NB;    // ScaLAPACK documentation
    const int LW   = (NB+BWU)*(BWL+BWU)+6*(BWL+BWU)*(BWL+2*BWU) + max(NRHS*(NB+2*BWL+4*BWU), 1); 
  
    double* A    = new double[LA];   // Matrix banded storage
    int*    ipiv = new int   [NB];   // Pivoting array
    double* x    = new double[NB];   // In: RHS vector, Out: Solution
    double* work = new double[LW];   // Workspace

    int desca[7];             // Descriptor for banded matrix
    desca[0] = 501;           // Type
    desca[1] = ctx;           // Context
    desca[2] = N;             // Problem size
    desca[3] = NB;            // Blocking of matrix
    desca[4] = 0;             // Process row/column
    desca[5] = 1+2*BWL+2*BWU; // Local leading dim
    desca[6] = 0;             // Reserved

    int descb[7];             // Descriptor for RHS
    descb[0] = 502;           // Type
    descb[1] = ctx;           // Context
    descb[2] = N;             // Problem size
    descb[3] = NB;            // Blocking of matrix
    descb[4] = 0;             // Process row/column
    descb[5] = NB;            // Local leading dim
    descb[6] = 0;             // Reserved

    // Populate banded A and x here 
    
    // ... Set up CBLACS grid

    // Perform the parallel solve.
    F77NAME(pdgbsv) (N, BWL, BWU, NRHS, A, JA, desca, ipiv, &x[0], IB, descb, work, LW, info);

    // Verify it completed successfully.
    if (info) {
    cout << "Error occurred in PDGBTRS: " << info << endl;
    }
    //---------------------------------------------------------------------------------------------------------

    // Finalize CBLACS and clean up memory
    //---------------------------------------------------------------------------------------------------------

    // pass

    //---------------------------------------------------------------------------------------------------------

    // Banded matrix solver

    //----------------------------------------------------------------------------------------------------------------
    // Parameters to be passed
    int n = (Nx - 2) * (Ny - 2);       // Problem size
    int kl = Nx - 2;                   // Lower diagonal bandwidth  (Determined solely by number of columns)
    int ku = kl;                       // Upper diagonal bandwidth (Ku = Kl always here)
    int nrhs = 1;                      // Number of RHS vectors (always 1 here)
    int ldab = 1 + 2*kl + ku;          // Number of rows in compressed matrix
    int ldb = n;                       // Size of RHS vector
    int info;                          // Additional possible info
    A = new double[ldab*n];            // Initialise banded array/matrix
    piv = new int[n];                  // Pivot data
    b = new double[n];                 // Output vector (vorticities)

    //----------------------------------------------------------------------------------------------------------------
    // Populate banded A matrix (column major format)
    for (int i=0; i<ldab*n; i++) {
        if ((i - 2*ku)%ldab == 0) {
            A[i] = 2/(dx*dx) + 2/(dy*dy);
        }
        else if ( ((i-ku)%ldab==0 && i > ku*ldab+1) || ((i+1)%ldab==0 && i < (n-ku)*ldab) ) {
            A[i] = -1/(dx*dx);
        }
        else if ( (i<ldab*n-ku && (i-2*ku-1)%ldab == 0 && (i+ku)%(ldab*ku)!=0) ||
                  (i>2*ku-1 && (i-2*ku+1)%ldab ==0 && (i-2*ku+1)%(ldab*ku)!=0) ) {
            A[i] = -1/(dy*dy);
        }
        else {
            A[i] = 0;
        }
    }

    // A Matrix Visualisation (displayed in column major format)
    // if (MPI_ROOT) {
    //     ofstream myfile6;
    //     myfile6.open("A_matrix.txt");
    //     for (int i=0; i<n; i++){
    //         for (int j=0; j<ldab; j++) {
    //             if (A[i*ldab + j] != 0) {
    //                 myfile6 << A[i*ldab + j] << " ";
    //             } 
    //             else {
    //                 myfile6 << A[i*ldab + j] << "   ";
    //             }
    //         }
    //         myfile6 << endl;
    //     }
    //     myfile6.close();
    // }


    //----------------------------------------------------------------------------------------------------------------
    // Populate b vector from passed in vorticity matrix argument
    // Convert the vorticity matrix into a long vector (exclude vorticities in the boundaries)
    // Incremement index by column, then by row
    double* vorticity_vec[n];
    // ofstream myfile5;
    // myfile5.open("b_vector.txt");
    for (int i=1; i<Ny-1; i++) {
        for (int j=1; j<Nx-1; j++) {
            vorticity_vec[(Ny*(i-1))+(j-1)] = (omega_new +i*Nx +j);
            b[(Ny-2)*(i-1)+(j-1)] = *vorticity_vec[(Ny*(i-1))+(j-1)];
            // myfile5 << b[(Ny-2)*(i-1)+(j-1)] << endl;
        }
    }
    // myfile5.close();

    // Visualise new internal streamfunction
    // if (MPI_ROOT) {
    //     ofstream myfile7;
    //     myfile7.open("stream_vector_new.txt");
    //     for (int i=0; i<(Ny-2)*(Nx-2); i++) {
    //         myfile7 << b[i] << endl;
    //     }
    //     myfile7.close();
    // }

}

void PoissonSolver::ReturnStream(double * psi_new, int Nx, int Ny) {
    cblas_dcopy((Nx-2)*(Ny-2), b, 1, psi_new, 1);
}