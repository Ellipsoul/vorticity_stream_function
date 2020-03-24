// Header Files
#include "LidDrivenCavity.h"
#include "PoissonSolver.h"

// Basic Libraries
#include <iostream>
#include <exception>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <typeinfo>
#include <chrono>

// External Libraries
#include <cblas.h>
#include <mpi.h>

using namespace std;

/// Class Constructor
LidDrivenCavity::LidDrivenCavity()
{
    // Get number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    np = world_size;

    // Get rank of process
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    mpiroot = (mpirank == 0);

    // Announce creation of instance
    if (mpiroot) {
        cout << "LidDrivenCavity instance created" << endl;
    }
}

/// Class Destructor
LidDrivenCavity::~LidDrivenCavity()
{
    delete[] psi_new;
}

/**
 * Verifies and sets the cavity domain
 * @param xlen     Horizontal dimension of cavity
 * @param ylen     Vertical dimension of cavity
 */ 
void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    // Catch and return errors if either value is not positive
    if (xlen <= 0) {
        throw domain_error("Lx must be a positive value. Program terminated");
    }
    if (ylen <= 0) {
        throw domain_error("Ly must be a positive value. Program terminated");
    }

    // Set values if checks pass
    Lx = xlen;
    Ly = ylen;
}

/**
 * Verifies and sets the discretisation grid size
 * @param nx        Grid points in x-direction
 * @param ny        Grid points in y-direction
 */ 
void LidDrivenCavity::SetGridSize(double nx, double ny)
{
    // Catch and return erros if either value is not a positive integer
    if (nx - int(nx) != 0 || nx <= 0) {
        throw domain_error("Nx must be a positive integer. Program terminated");
    }
    if (ny - int(ny) != 0 || ny <= 0) {
        throw domain_error("Ny must be a positive integer. Program terminated");
    }

    // Set values if checks pass
    Nx = int(nx);
    Ny = int(ny);
}

/**
 * Verifies and sets the simulation duration
 * @param finalt        Simulation duration
 */ 
void LidDrivenCavity::SetFinalTime(double finalt)
{
    // Catch and return error if time input is not positive
    if (finalt <= 0) {
        throw domain_error("T must be a positive value. Program terminated");
    }

    // Set value if check passes
    T = finalt;
}

/**
 * Verifies and sets the Reynolds Number
 * @param re        Reynolds number
 */ 
void LidDrivenCavity::SetReynoldsNumber(double re)
{
    // Catch and return error if Reynolds number is not positive
    if (re <= 0) {
        throw domain_error("Re must be a positive value. Program terminated");
    }

    // Set value if check passes
    Re = re;
}

/**
 * Verifies and sets the time increment step
 * @param deltat        Time increment step
 * @param Lx            Horizontal dimension of cavity
 * @param Ly            Vertical dimension of cavity
 * @param Nx            Grid points in x-direction
 * @param Ny            Grid points in y-direction
 * @param Re            Reynolds number
 */ 
void LidDrivenCavity::SetTimeStep(double deltat, double Lx, double Ly, unsigned int Nx, unsigned int Ny, double Re)
{
    // Catch and return error if timestep is not positive
    if (deltat <= 0 || deltat >= Re*(Lx/Nx)*(Ly/Ny)/4) {
        throw domain_error("dt must be a positive value and within the allowed range. Program terminated");
    }

    // Set value if check passes
    dt = deltat;
}

/**
 * Verifies and sets the processor partitions for parallel programming
 * @param px        Number of partitions in x-direction
 * @param py        Number of partitions in y-direction
 * @param nx        Grid points in x-direction
 * @param ny        Grid points in y-direction
 */ 
void LidDrivenCavity::SetPartitions(double px, double py, int nx, int ny)
{
    // Catch and return errors if either value is not a positive integer, or fails to properly divide domain
    if (px - int(px) != 0 || px <= 0 || (nx - 1)/px != int((nx - 1)/px)) {
        throw domain_error("Px must be a positive integer and appropriately divide the domain. Program terminated");
    }
    if (py - int(py) != 0 || py <= 0 || (ny - 1)/py != int((ny - 1)/py)) {
        throw domain_error("Py must be a positive integer and appropriately divide the domain. Program terminated");
    }

    // Set values if checks pass
    Px = int(px);
    Py = int(py);
}

/**
 * Verifies the number of processors is consistent with partitions
 * @param px        Number of partitions in x-direction
 * @param py        Number of partitions in y-direction
 */ 
void LidDrivenCavity::VerifyProcessors(double px, double py) {
    // Catch and return error if there are too many processors for the current domain size
    if (np > (Nx-2)*(Ny-2)/(2*Nx-3)) {
        throw domain_error("Too many processors entered for an efficient parallel solution. Program terminated");
    }
    // Catch and return error if number of processors entered fails to fill partition grid
    if (np != px*py) {
        throw domain_error("Number of processors provided must exactly fill the partition grid. Program terminated");
    }
}

/**
 * Performs all the above class functions and prepares the solver
 * @param Lx_arg        Horizontal dimension of cavity
 * @param Ly_arg        Vertical dimension of cavity
 * @param Nx_arg        Grid points in x-direction
 * @param Ny_arg        Grid points in y-direction
 * @param Px_arg        Number of partitions in x-direction
 * @param Py_arg        Number of partitions in y-direction
 * @param dt_arg        Time increment step
 * @param T_arg         Simulation duration
 * @param Re_arg        Reynolds number
 */
void LidDrivenCavity::Verify(const double Lx_arg, const double Ly_arg, const double Nx_arg, const double Ny_arg,
                                const double Px_arg, const double Py_arg, const double dt_arg, const double T_arg, 
                                const double Re_arg) 
{
    // Set cavity/domain size
    try {
        SetDomainSize(Lx_arg, Ly_arg);
    }
    catch (domain_error& e) {
        if (mpiroot) {
            cout << "Error: " << e.what() << endl;
        }
        exit(EXIT_FAILURE);
    }

    // Set grid size
    try {
        SetGridSize(Nx_arg, Ny_arg);
    }
    catch (domain_error& e) {
        if (mpiroot) {
            cout << "Error: " << e.what() << endl;
        }
        exit(EXIT_FAILURE);
    }

    // Set time step
    try {
        SetTimeStep(dt_arg, Lx_arg, Ly_arg, Nx_arg, Ny_arg, Re_arg);
    }
    catch (domain_error& e) {
        if (mpiroot) {
            cout << "Error: " << e.what() << endl;
        }
        exit(EXIT_FAILURE);
    }

    // Set final time
    try {
        SetFinalTime(T_arg);
    }
    catch (domain_error& e) {
        if (mpiroot) {
            cout << "Error: " << e.what() << endl;
        }
        exit(EXIT_FAILURE);
    }

    // Set Reynolds number
    try {
        SetReynoldsNumber(Re_arg);
    }
    catch (domain_error& e) {
        if (mpiroot) {
        cout << "Error: " << e.what() << endl;
        }
        exit(EXIT_FAILURE);
    }

    // Set number of partitions
    try {
        SetPartitions(Px_arg, Py_arg, Nx_arg, Ny_arg);
    }
    catch (domain_error& e) {
        if (mpiroot) {
            cout << "Error: " << e.what() << endl;
        }
        exit(EXIT_FAILURE);
    }

    // Verify number of processors
    try {
        VerifyProcessors(Px_arg, Py_arg);
    }
    catch (domain_error& e) {
        if (mpiroot) {
            cout << "Error: " << e.what() << endl;
        }
        exit(EXIT_FAILURE);
    }
}

/**
 * Executes the solver
 */
void LidDrivenCavity::Solve()
{  
    // Start runtime timer
    auto start = chrono::high_resolution_clock::now();

    // Initialise initial condition matrices
    const int Nx_const = Nx;
    const int Ny_const = Ny;
    double psi[Nx_const][Ny_const];
    double omega[Nx_const][Ny_const];
    double omega_new[Nx_const][Ny_const];

    // Initial vorticity and stream-function zero matrices
    for (int i=0; i<Nx; i++){
        for (int j=0; j<Ny; j++) {
            psi[i][j] = 0;
            omega[i][j] = 0;
            omega_new[i][j] = 0;
        }
    }

    // Defining helper variables and lid velocity
    double dx = Lx/(Nx-1);
    double dy = Ly/(Ny-1);
    double U = 1.0;
    int t_steps = ceil(T/dt);

    // Create new instance of the Poisson Solver
    PoissonSolver* poisson = new PoissonSolver();

    // Loop through each time increment
    for (int i=1; i<t_steps; i++) { 
        
        // Calculating vorticity boundary conditions at time t
        //---------------------------------------------------------------------------------------------------------
        for (int j=0; j<Nx; j++) {
            omega[0][j] = 2/(dy*dy) * (psi[0][j] - psi[1][j]) - 2*U/dy;  // Top Surface
            omega_new[0][j] = omega[0][j];
            omega[Nx-1][j] = 2/(dy*dy) * (psi[Nx-1][j] - psi[Nx-2][j]);  // Bottom Surface
            omega_new[Nx-1][j] = omega[Nx-1][j];
        }

        for (int j=1; j<Ny; j++) {
            omega[j][0] = 2/(dx*dx) * (psi[j][0] - psi[j][1]);           // Left Surface
            omega_new[j][0] = omega[j][0];
            omega[j][Nx-1] = 2/(dx*dx) * (psi[j][Nx-1] - psi[j][Nx-2]);  // Right Surface
            omega_new[j][Nx-1] = omega[j][Nx-1];
        }
        //---------------------------------------------------------------------------------------------------------

        // Calculate interior vorticity at time t
        //---------------------------------------------------------------------------------------------------------
        for (int j=1; j<Nx-1; j++) {
            for (int k=1; k<Ny-1; k++) {
                omega[j][k] = - ( (psi[j-1][k] - 2*psi[j][k] + psi[j+1][k])/(dx*dx) + 
                                  (psi[j][k+1] - 2*psi[j][k] + psi[j][k-1])/(dy*dy) );
            }
        }
        //---------------------------------------------------------------------------------------------------------

        // Calculate interior vorticity at time t + dt
        //---------------------------------------------------------------------------------------------------------
        for (int j=1; j<Nx-1; j++) {
            for (int k=1; k<Ny-1; k++) {
                omega_new[j][k] = ( (1/Re)* ( (omega[j-1][k]-2*omega[j][k]+omega[j+1][k])/(dx*dx) + 
                                    (omega[j][k+1]-2*omega[j][k]+omega[j][k-1])/(dy*dy) ) + 
                                    (psi[j-1][k]-psi[j+1][k])/(2*dx) * (omega[j][k+1]-omega[j][k-1])/(2*dy) -
                                    (psi[j][k+1]-psi[j][k-1])/(2*dy) * (omega[j-1][k]-omega[j+1][k])/(2*dx) ) * dt +
                                    omega[j][k];
            }
        }
        //---------------------------------------------------------------------------------------------------------

        // Vorticity and Stream-function visualisation (uncomment if needed)
        // if (mpiroot) {
        //     ofstream myfile1;
        //     ofstream myfile2;
        //     ofstream myfile3;
        //     myfile1.open("stream_matrix_old.txt");
        //     myfile2.open("vorticity_matrix_old.txt");
        //     myfile3.open("vorticity_matrix_new.txt");
        //     for (int i=0; i<Nx; i++){
        //         for (int j=0; j<Ny; j++) {
        //             myfile1 << psi[i][j] << " ";
        //             myfile2 << omega[i][j] << " ";
        //             myfile3 << omega_new[i][j] << " ";
        //         }
        //         myfile1 << endl;
        //         myfile2 << endl;
        //         myfile3 << endl;
        //     }
        //     myfile1.close();
        //     myfile2.close();
        //     myfile3.close();
        // }

        // Solve the Poisson problem to calculate stream-function at time t + dt
        //---------------------------------------------------------------------------------------------------------        
        // Execute poisson solver
        poisson -> SolvePoisson((double*)omega_new, Ny, Nx, dx, dy);

        // Retrieve values for the new stream-function
        psi_new = new double[(Nx-2)*(Ny-2)];
        poisson -> ReturnStream(psi_new, Nx, Ny);
        //---------------------------------------------------------------------------------------------------------

        // Update the the stream-function into the original psi matrix
        //--------------------------------------------------------------------------------------------------------- 
        for (int i=1; i<Ny-1; i++) {
            for (int j=1; j<Nx-1; j++) {
                psi[i][j] = psi_new[(Ny-2)*(i-1) + (j-1)];
            }
        }
        //---------------------------------------------------------------------------------------------------------

        // Updated stream-function matrix visualisation (uncomment if needed)
        // if (mpiroot) {
        //     ofstream myfile8;
        //     myfile8.open("stream_matrix_new.txt");
        //     for (int i=0; i<Nx; i++) {
        //         for (int j=0; j<Ny; j++) {
        //             myfile8 << psi[i][j] << " ";
        //         }
        //         myfile8 << endl;
        //     }
        //     myfile8.close();
        // }
    }

    // Final streamfunction matrix found, solving for velocities and output to a file
    //-------------------------------------------------------------------------------------------------------------
    if (mpiroot) {
        double u[Ny-1][Nx];
        double v[Nx][Ny-1];

        // u velocity
        ofstream u_file;
        u_file.open("u_velocity.txt");
        for (int i=0; i<Ny-1; i++) {
            for (int j=0; j<Nx; j++) {
                u[i][j] = (psi[i][j] - psi[i+1][j])/dy;
                u_file << u[i][j] << " ";
            }
            u_file << endl;
        }
        u_file.close();

        // v velocity
        ofstream v_file;
        v_file.open("v_velocity.txt");
        for (int i=0; i<Ny; i++) {
            for (int j=0; j<Nx-1; j++) {
                v[i][j] = (psi[i][j+1] - psi[i][j])/dx;
                v_file << v[i][j] << " ";
            }
            v_file << endl;
        }
        v_file.close();

        // Final streamfunction and vorticity matrices
        ofstream psi_file;
        ofstream omega_file;
        psi_file.open("psi_final.txt");
        omega_file.open("omega_final.txt");
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny; j++) {
                psi_file << psi[i][j] << " ";
                omega_file << omega[i][j] << " ";
            }
            psi_file << endl;
            omega_file << endl;
        }
        psi_file.close();
        omega_file.close();
    }
    //-------------------------------------------------------------------------------------------------------------
    // Stop timer and calculate code runtime
    auto stop = chrono::high_resolution_clock::now();
    if (mpiroot) {
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Elapsed execution duration: " << duration.count() << " milliseconds" << endl;
    }

    // Run destructor for poisson class instance
    poisson -> ~PoissonSolver();
}
