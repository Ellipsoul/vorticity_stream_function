# Define required version of CMake (check current version)
cmake_minimum_required(VERSION 3.16)

# Name of project (Keep same, will be the name of output file)
project(lidDrivenCavity)

# Setting to C++11 version
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Have the project name, then all source files following that
# Header files will automatically be included with the linked cpp file
add_executable(
    lidDrivenCavity
    LidDrivenCavity.cpp
    LidDrivenCavitySolver.cpp
    PoissonSolver.cpp
)

# Attach all libraries here
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(SCALAPACK REQUIRED)
find_package(MPI REQUIRED)

target_include_directories(lidDrivenCavity PRIVATE ${BLAS_INCLUDE_DIR})
target_include_directories(lidDrivenCavity PRIVATE ${Boost_INCLUDE_DIR})
target_include_directories(lidDrivenCavity PRIVATE ${LAPACK_INCLUDE_DIR})
target_include_directories(lidDrivenCavity PRIVATE ${MPI_INCLUDE_PATH})

# Link the corresponding libraries from the 'target' to the found library
# target_link_libraries(myprog ${BLAS_LIBRARIES})
target_link_libraries(lidDrivenCavity ${BLAS_LIBRARIES} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_LIBRARIES})

# Instructions for creating and running CMake files:
#   1. Filename is "CMakeLists.txt"
#   2. Build in a separate directory:
#       mkdir build
#       cd build
#       cmake ../
#   3. In this new directory run the Makefile:
#       make
#   4. The program executable will be made with the project name. Run it like an ./a.out file (without extension)
#       ./lidDrivenCavity

# Compile and execute command from project directory:
# mkdir build; cd build; cmake ../; make; ./lidDrivenCavity --Lx 1.0 --Ly 1.0 --Nx 5 --Ny 5 --Px 4 --Py 4 --dt 0.01 --T 10 --Re 4000

# Full command line recompile and execute:
# make clean;cd ..;rmdir build;mkdir build;cd build;cmake ../;make;./lidDrivenCavity --Lx 1.0 --Ly 1.0 --Nx 5 --Ny 5 --Px 4 --Py 4 --dt 0.01 --T 10 --Re 4000

# Pass in these arguments when running
# ./lidDrivenCavity --Lx 1.0 --Ly 1.0 --Nx 5 --Ny 5 --Px 4 --Py 4 --dt 0.05 --T 10 --Re 400

# Compile in command line if packages fail
# g++ LidDrivenCavity.cpp LidDrivenCavitySolver.cpp -lblas -lboost_program_options -std=c++11