### Vorticity Stream Function

Parallel numeric code for solving the vorticity stream function of incompressible Navier-Stokes Equations. Written as part of a High Performance Computing assignment for 3rd Year Aeronautical Engineering at Imperial College London.

##### Running the solver

To run the solver in the command line from the main directory:

``` C++
mkdir build; cd build; cmake ../; make; ./lidDrivenCavity --Lx 1.0 --Ly 1.0 --Nx 5 --Ny 5 --Px 4 --Py 4 --dt 0.01 --T 10 --Re 4000
```

