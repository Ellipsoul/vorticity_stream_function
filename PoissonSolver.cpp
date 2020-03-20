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

// Setting up BLACS grid
extern "C" {
    // CBlacs Declarations
    void Cblacs_pinfo(int*, int*);
    void Cblacs_get(int, int, int*);
    void Cblacs_gridinit(int*, const char*, int, int);
    void Cblacs_pcoord(int, int, int*, int*);
    void Cblacs_gridexit(int);
    void Cblacs_barrier(int, const char*);
    void Cdgerv2d(int, int, int, double*, int, int, int);
    void Cdgesd2d(int, int, int, double*, int, int, int);
 
    int numroc_(int*, int*, int*, int*, int*);
}

// Defining LAPACK routine for solving banded linear system
#define F77NAME(x) x##_
extern "C" {
    void F77NAME(dgbsv) (const int& n, const int& kl, const int& ku, const int& nrhs, const double * A,
                         const int& ldab, int * piv, double * B, const int &ldb, int& info);
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

    //----------------------------------------------------------------------------------------------------------------
    // Running the solver
    F77NAME(dgbsv) (n, kl, ku, nrhs, A, ldab, piv, b, ldb, info);

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