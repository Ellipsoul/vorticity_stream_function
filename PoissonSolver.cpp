// Header Files
#include "PoissonSolver.h"

// Basic Libraries
#include <iostream>
#include <math.h>
#include <fstream>
#include <cstdlib>

// External Libraries
#include <cblas.h>
#include <mpi.h>

using namespace std;

// Declaring Blacs function
extern "C" {
    // Output parameters have a * while input parameters do not
    void Cblacs_pinfo(int*, int*);
    void Cblacs_get(int, int, int*);
    void Cblacs_gridinit(int*, const char*, int, int);
    void Cblacs_gridinfo(int, int*, int*, int*, int*);
    void Cblacs_barrier(int, const char*);
    void Cblacs_gridexit(int);
    void Cblacs_exit(int);
}

// Defining LAPACK routine for solving parallel banded linear system
#define F77NAME(x) x##_
extern "C" {
    void F77NAME(pdgbsv) (const int& N, const int& BWL, const int& BWU, const int& NRHS, double * A,
                          const int& JA, int * desca, int * ipiv, double * x, const int& IB, int * descb,
                          double * work, const int& LW, int& info) ;
}

/// Class Constructor
PoissonSolver::PoissonSolver()
{
    // Declare MPI variables
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    mpiroot = (mpirank == 0);
    if (mpiroot) {
        cout << "Poisson Solver instance created" << endl;
    }

    // Initialise CBLACS
    //----------------------------------------------------------------------------------------------------------------
    // Initialises the BLACS world communicator
    // npe: total number of processes
    // mype: process rank (starting from 0)
    Cblacs_pinfo(&mype, &npe);

    // Get the default system context (i.e. MPI_COMM_WORLD)
    Cblacs_get( 0, 0, &ctx );

    // Initialise a process grid of 1 rows and npe columns
    Cblacs_gridinit( &ctx, &order, 1, npe );

    // Get info about the grid to verify it is set up
    // myrow and mycol are indexed from 0
    Cblacs_gridinfo( ctx, &nrow, &ncol, &myrow, &mycol);
    //----------------------------------------------------------------------------------------------------------------
}

/// Class Destructor
PoissonSolver::~PoissonSolver()
{
    delete[] A;
    delete[] b;

    delete[] A_loc;
    delete[] x;
    delete[] ipiv;
    delete[] work;

    Cblacs_gridexit(ctx);
    Cblacs_exit(0);
}

/**
 * Executes the Poisson solver
 * @param omega_new     Vorticity matrix
 * @param Ny            Grid points in y-direction
 * @param Nx            Grid points in x-direction
 * @param dx            Horizontal spacial increment
 * @param dy            Vertical spacial increment
 */ 
void PoissonSolver::SolvePoisson(double* omega_new, int Ny, int Nx, double dx, double dy) {
    
    // Visualise passed vorticity matrix (uncomment as needed)
    // if (mpiroot) {
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

    // Populate Global A matrix (column major format)
    //----------------------------------------------------------------------------------------------------------------
    n = (Nx-2)*(Ny-2);             // Number of columns of global A matrix
    kl = Nx - 2;                   // Lower diagonal bandwidth  (Determined solely by number of columns)
    ku = kl;                       // Upper diagonal bandwidth (Ku = Kl always here)
    ldab = 1 + 2*kl + 2*ku;        // Number of rows in banded global matrix

    A = new double[ldab*n];        // Initialise banded array/matrix

    // Populate Global A matrix
    for (int i=0; i<ldab*n; i++) {
        if ((i - 3*ku)%ldab == 0) {
            A[i] = 2/(dx*dx) + 2/(dy*dy);
        }
        else if ( ((i-2*ku)%ldab==0 && i > ku*ldab+1) || ((i+1)%ldab==0 && i < (n-ku)*ldab) ) {
            A[i] = -1/(dx*dx);
        }
        else if ( (i<ldab*n-ku && (i-3*ku-1)%ldab == 0 && (i+ku)%(ldab*ku)!=0) ||
                (i>3*ku-1 && (i-3*ku+1)%ldab==0 && (i-3*ku+1)%(ldab*ku)!=0) ) {
            A[i] = -1/(dy*dy);
        }
        else {
            A[i] = 0;
        }
    }
    //----------------------------------------------------------------------------------------------------------------

    // A Matrix Visualisation (displayed in column major format, uncomment if needed)
    // if (mpiroot) {
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
    
    // Populate global b vector from passed in vorticity matrix argument
    //----------------------------------------------------------------------------------------------------------------

    // Convert the vorticity matrix into a long vector (exclude vorticities in the boundaries)
    // Incremement index by column, then by row
    double* vorticity_vec[n];
    b = new double[n];

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

    // Parallel Banded Matrix Solver
    //----------------------------------------------------------------------------------------------------------------
    int info;                                   // Status value
    const int N    = (Nx-2)*(Ny-2);             // Total problem size
    const int NB   = ceil(1.0*N/npe);           // Blocking size (number of columns per process)
    const int BWL  = Nx-2;                      // Lower bandwidth
    const int BWU  = Nx-2;                      // Upper bandwidth
    const int NRHS = 1;                         // Number of RHS to solve
    const int JA   = 1;                         // Start offset in matrix (to use just a submatrix)
    const int IB   = 1;                         // Start offset in RHS vector (to use just a subvector)
    const int LA   = (1 + 2*BWL + 2*BWU)*NB;    // ScaLAPACK documentation
    const int LW   = (NB+BWU)*(BWL+BWU)+6*(BWL+BWU)*(BWL+2*BWU) + max(NRHS*(NB+2*BWL+4*BWU), 1); 

    A_loc = new double[LA];   // Matrix banded storage
    ipiv  = new int   [NB];   // Pivoting array
    x     = new double[NB];   // In: RHS vector, Out: Solution
    work  = new double[LW];   // Workspace

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

    // Local A matrix
    //----------------------------------------------------------------------------------------------------------------
    // Populate local A matrix
    for (int i=0; i<LA; i++) {
        if ((mype+1) != npe && N%NB != 0) {  // All processes except last (if padding required)
            A_loc[i] = A[i + mype*(1 + 2*BWL + 2*BWU)*NB];
        }
        else {
            if (i < N*(1 + 2*BWL + 2*BWU)) {
                A_loc[i] = A[i + mype*(1 + 2*BWL + 2*BWU)*NB];
            }
            else {
                A_loc[i] = 0;
            }
        }
    }
    //----------------------------------------------------------------------------------------------------------------

    // Local A matrix visualisation (uncomment as needed)
    // if (mpiroot) {
    //     ofstream A_locfile;
    //     A_locfile.open("Local_A.txt");
    //     for (int i=0; i<NB; i++) {
    //         for (int j=0; j<(1+2*BWL+2*BWU); j++) {
    //             A_locfile  << A_loc[(1+2*BWL+2*BWU)*i + j] << " ";
    //         }
    //         A_locfile << endl;
    //     }
    //     A_locfile.close();
    // }

    // Local x vector
    //----------------------------------------------------------------------------------------------------------------
    // Populate local x vector
    for (int i=0; i<NB; i++) {
        if ((mype+1) != npe && N%NB != 0) {  // All processes except last (if padding required)
            x[i] = b[i + mype*(NB)];
        }
        else {
            if (i < N) {
                x[i] = b[i + mype*(NB)];
            }
            else {
                x[i] = 0;
            }
        }
    }
    //----------------------------------------------------------------------------------------------------------------

    // Local x vector visualisation (uncomment as needed)
    // if (!mpiroot) {
    //     ofstream x_locfile;
    //     x_locfile.open("Local_x.txt");
    //     for (int i=0; i<NB; i++) {
    //         x_locfile << x[i] << endl;
    //     }
    //     x_locfile.close();
    // }
    
    // Perform parallel solve
    //---------------------------------------------------------------------------------------------------------
    F77NAME(pdgbsv) (N, BWL, BWU, NRHS, A, JA, desca, ipiv, &x[0], IB, descb, work, LW, info);
    // Verify it completed successfully.
    if (info) {
    cout << "Error occurred in PDGBTRS: " << info << endl;
    }
    //---------------------------------------------------------------------------------------------------------
    
    // Visualise new internal streamfunction (uncomment as needed)
    // if (mpiroot) {
    //     ofstream myfile7;
    //     myfile7.open("stream_vector_new.txt");
    //     for (int i=0; i<(Ny-2)*(Nx-2); i++) {
    //         myfile7 << b[i] << endl;
    //     }
    //     myfile7.close();
    // }

}

/**
 * Returns the calculated stream-function
 * @param psi_new       Stream-function vector
 * @param Nx            Grid points in x-direction
 * @param Ny            Grid points in y-direction
 */ 
void PoissonSolver::ReturnStream(double * psi_new, int Nx, int Ny) {
    // Assemble a global vorticity vector
    MPI_Allgather(x, (Nx-2)*(Ny-2), MPI_DOUBLE, b, (Nx-2)*(Ny-2), MPI_DOUBLE, MPI_COMM_WORLD);

    // Copy global streamfunction vector back to LidDrivenCavity
    cblas_dcopy((Nx-2)*(Ny-2), b, 1, psi_new, 1);
}