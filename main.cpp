// Build command:
// g++ main.cpp -o main
// g++ is the compiler, main.cpp is the input file, and "-o main" specifies that the executable file is called "main".
// TODO: When the connectors reconnect, they "bounce" to a upright position
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <vector>
#include "headers/System.h"
#include "headers/Block.h"
#include "headers/Params.h"

using namespace std;


// Forward declare functions
void writeVectorToFile(ofstream & outFile, vector<double> &vec, int numBlocks);

int main() // This function runs when you execute the program.
{
	// They constructor takes care of the parameters
    Params params("params.txt");
	System system(params);
	const int writeFrequency = 1000 * system.m_tStop;
	cout << "The write frequency is " << writeFrequency << endl;

	// Create output streams. These are closed upon the deletion of System
	system.openPositionsFile("output/positions.bin");
	system.openStatesFile("output/states.bin");
	system.openForcesFile("output/forces.bin");
	system.openConnectorsFile("output/connectorForces.bin");

	clock_t start, end;
	start = clock();

	int counter = 0;
    system.m_t = 0;
	while (system.m_t<system.m_tStop)
	{
        // Simulate a timestep of the entire system
		system.simulate();

		// modulo operation to check whether to write output to file on this timestep
		if ( (counter%writeFrequency) == 0)
            system.dumpData();

		system.m_t += system.m_dt;
		counter ++;
	}
    end = clock();

	// Output parameters to file
    ofstream outFileParameters("output/parameters.txt");
	outFileParameters << system;
	outFileParameters.close();

	cout << "Ran " << counter << " integration steps for "<<system.m_numBlocksY
         <<"x" << system.m_numBlocksX << " blocks with " << system.m_numConnectors << " micro-junctions each in "
         << ((double)end - (double)start)/CLOCKS_PER_SEC << " seconds" << endl;
	return 0;

}


void writeVectorToFile(ofstream & outFile, vector<double> &vec, int numBlocks)
{
	outFile.write(reinterpret_cast<char*>(&vec[0]), numBlocks*sizeof(double));
}

