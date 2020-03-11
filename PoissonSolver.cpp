// Header Files
#include "PoissonSolver.h"

// Basic
#include <iostream>
#include <math.h>
#include <fstream>

// External Libraries
#include "cblas.h"

using namespace std;

// Constructor
PoissonSolver::PoissonSolver()
{
}

// Destructor
PoissonSolver::~PoissonSolver()
{
}

void PoissonSolver::SolvePoisson(double * omega_new, int Ny, int Nx) {
    ofstream myfile4;
    myfile4.open("omega_trans.txt");
    for (int i=0; i<Ny; i++) {
        for (int j=0; j<Nx; j++) {
            myfile4 << *(omega_new + i*Nx + j) << " ";
        }
        myfile4 << endl;
    }
    myfile4.close();
}