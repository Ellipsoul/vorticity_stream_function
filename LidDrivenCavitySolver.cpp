#include <iostream>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    // Create a new instance of the LidDrivenCavity class
    // Uses a pointer so solver is actually storing the address of the new instance
    LidDrivenCavity* solver = new LidDrivenCavity();
    
    // Configure the solver here...
    // Arrow dereferences and initialises the solver at that address
    solver -> Initialise();

    // Run the solver
    solver -> Integrate();

	return 0;
}