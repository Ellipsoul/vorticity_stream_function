### Vorticity Stream Function

Parallel numeric code for solving the vorticity stream function of incompressible Navier-Stokes Equations. Written as part of a High Performance Computing assignment for 3rd Year Aeronautical Engineering at Imperial College London.

GitHub Repository: https://github.com/Ellipsoul/vorticity_stream_function
Branch "Submission_SpitFire" was submitted for this assignment

##### Running the solver

To run the solver in the command line from the main directory:

``` C++
mkdir build; cd build; cmake ../; make; mpiexec -np 1 ./lidDrivenCavity --Lx 1.0 --Ly 1.0 --Nx 15 --Ny 15 --Px 1 --Py 1 --dt 0.0005 --T 5 --Re 100
```

Note: As I was unable to fully implement the parallel solver (further explanation within the code), please set parameters Px, Py and np to 1 so that the code will be able to run in serial.
