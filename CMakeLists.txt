# Define required version of CMake (check current version)
cmake_minimum_required(VERSION 3.16)

# Name of project (Keep same, will be the name of output file)
project(lidDrivenCavity)

# Setting to C++11 version
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories("${CMAKE_SOURCE_DIR}/boost_1_72_0")

# Have the project name, then all source files following that
# Header files will automatically be included with the linked cpp file
add_executable(
    lidDrivenCavity
    LidDrivenCavity.cpp
    LidDrivenCavitySolver.cpp
)

# Attach all libraries here
find_package(Boost REQUIRED)

# Link the corresponding libraries from the 'target' to the found library
# target_link_libraries(myprog ${BLAS_LIBRARIES})
target_link_libraries(lidDrivenCavity Boost)

# Instructions for creating and running CMake files:
#   1. Filename is "CMakeLists.txt"
#   2. Build in a separate directory:
#       mkdir build
#       cd build
#       cmake ../
#   3. In this new directory run the Makefile:
#       make
#   4. The program executable will be made with the project name. Run it like an ./a.out file (without extension)
#       ./myprog