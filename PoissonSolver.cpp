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
}

// Function to solve the Poisson problem
void PoissonSolver::SolvePoisson(double * omega_new, int Ny, int Nx) {
    
    // Writes the contents of the passed vorticity matrix into a file
    ofstream myfile4;
    myfile4.open("vorticity_trans.txt");
    for (int i=0; i<Ny; i++) {
        for (int j=0; j<Nx; j++) {
            myfile4 << *(omega_new + i*Nx + j) << " ";
        }
        myfile4 << endl;
    }
    myfile4.close();

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
    double* A = new double[ldab*n];    // Initialise banded array/matrix
    int* piv = new int[n];             // Pivot data
    double* b = new double[n];         // Output vector (vorticities)

    //----------------------------------------------------------------------------------------------------------------
    // Populate banded A matrix
    for (int i=0; i<ldab*n; i++) {
        if ((i - 2*ku)%ldab == 0) {
            A[i] = 4;
        }
        else if ( ((i-ku)%ldab==0 & i > ku*ldab+1) || ((i+1)%ldab==0 & i < (n-ku)*ldab) || 
                  (i<ldab*n-ku & (i-2*ku-1)%ldab == 0 & (i+ku)%(ldab*ku)!=0) ||
                  (i>2*ku-1 & (i-2*ku+1)%ldab ==0 & (i-2*ku+1)%(ldab*ku)!=0) ) {
            A[i] = -1;
        }
        else {
            A[i] = 0;
        }
    }

    // A Matrix Visualisation
    ofstream myfile6;
    myfile6.open("A_matrix.txt");
    for (int i=0; i<n; i++){
        for (int j=0; j<ldab; j++) {
            if (A[i*ldab + j] == -1) {
                myfile6 << A[i*ldab + j] << " ";
            } 
            else {
                myfile6 << A[i*ldab + j] << "  ";
            }
        }
        myfile6 << endl;
    }
    myfile6.close();

    //----------------------------------------------------------------------------------------------------------------
    // Populate b vector from passed in vorticity matrix argument
    // Convert the vorticity matrix into a long vector (exclude vorticities in the boundaries)
    // Incremement index by column, then by row
    double* vorticity_vec[n];
    ofstream myfile5;
    myfile5.open("vorticity_vector.txt");
    for (int i=1; i<Ny-1; i++) {
        for (int j=1; j<Nx-1; j++) {
            vorticity_vec[(Ny*(i-1))+(j-1)] = (omega_new +i*Nx +j);
            b[(Ny-2)*(i-1)+(j-1)] = *vorticity_vec[(Ny*(i-1))+(j-1)];
            myfile5 << b[(Ny-2)*(i-1)+(j-1)] << endl;
        }
    }
    myfile5.close();

    //----------------------------------------------------------------------------------------------------------------
    // Running the solver
    F77NAME(dgbsv) (n, kl, ku, nrhs, A, ldab, piv, b, ldb, info);

    cout << b[1] << endl;

}

void PoissonSolver::PassPoisson(double* psi_new[]) {
    for (int i=0; i<n; i++) {
        psi_new[i] = &b[i];
    }
} 