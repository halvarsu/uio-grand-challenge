// Build command:
// g++ main.cpp -o main
// g++ is the compiler, main.cpp is the input file, and "-o main" specifies that the executable file is called "main".

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <vector>
#include "headers/Blocks.h"
#include "headers/Params.h"

using namespace std;


// Forward declare functions
void writeArrayToFile(ofstream & outFile, double * array, int numBlocks);
void writeVectorToFile(ofstream & outFile, vector<double> &vec, int numBlocks);

int main() // This function runs when you execute the program.
{
	// The constructor takes care of the parameters
    Params params("params.txt");
	Blocks blocks(params);
    cout << blocks;

	int writeFrequency		= (int) (1000 * blocks.m_tStop);
	cout << "The write frequency is " << writeFrequency << endl;

	// Create output streams
	ofstream outFilePositions("output/positions.bin");
	ofstream outFileStates("output/states.bin");
	ofstream outFileParameters("output/parameters.txt");
	ofstream outFileForces("output/forces.bin");
	ofstream outFileConnectorForces("output/connectorForces.bin");

	clock_t start, end;
	start = clock();

	int counter = 0;
	while (blocks.m_t<blocks.m_tStop)
	{
		// Calculate forces
		blocks.calculateForces();
		blocks.integrate();

		// modulo operation to check whether to write output to file on this timestep
		if ( (counter%writeFrequency) == 0)
		{
			writeArrayToFile(outFilePositions, blocks.getPositions(), blocks.m_numBlocks);
            writeArrayToFile(outFileStates, blocks.getStates(), blocks.m_numBlocks);
			writeArrayToFile(outFileForces, blocks.getForces(), blocks.m_numBlocks);
			writeArrayToFile(outFileConnectorForces, blocks.getConnectorForces(), blocks.m_numBlocks);
		}
		blocks.m_t += blocks.m_dt;
		counter ++;
	}
    end = clock();

	// Output parameters to file
	outFileParameters << blocks;

	// Close output files
	outFilePositions.close();
	outFileParameters.close();
	outFileStates.close();
	outFileForces.close();
	outFileConnectorForces.close();

	cout << "Ran " << counter << " integration steps for "<<blocks.m_numBlocks<<" blocks in " << ((double)end - (double)start)/
		CLOCKS_PER_SEC << " seconds" << endl;
	return 0;

}


void writeVectorToFile(ofstream & outFile, vector<double> &vec, int numBlocks)
{
	outFile.write(reinterpret_cast<char*>(&vec[0]), numBlocks*sizeof(double));
}

void writeArrayToFile(ofstream & outFile, double * array, int numBlocks)
{
    if(array)
        outFile.write(reinterpret_cast<char*>(array), numBlocks*sizeof(array[0]));
}
