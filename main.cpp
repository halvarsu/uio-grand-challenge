// Build command:
// g++ main.cpp -o main
// g++ is the compiler, main.cpp is the input file, and "-o main" specifies that the executable file is called "main".

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <vector>
#include "headers/Block.h"


using namespace std;


// Forward declare functions
void writeArrayToFile(ofstream & outFile, double * array, int numBlocks);
void writeVectorToFile(ofstream & outFile, vector<double> &vec, int numBlocks);



int main() // This function runs when you execute the program.
{
	// Choose parameters

	const int numBlocks		= 70;
	const double dt	   		= 1e-7;
	const double tStop 		= 0.01;

	// The constructor takes care of the parameters
	Block blocks(numBlocks);

	int writeFrequency		= 10;

	// Create output streams
	ofstream outFilePositions("output/positions.bin");
	ofstream outFileStates("output/states.bin");
	ofstream outFileParameters("output/parameters.txt");

	clock_t start, end;
	start = clock();

	int counter = 0;
	while (blocks.m_t<tStop)
	{
		// Calculate forces
		blocks.calculateForces(dt);
		blocks.integrate(dt);

		// modulo operation to check whether to write output to file on this timestep
		if ( (counter%writeFrequency) == 0)
		{
			writeArrayToFile(outFilePositions, blocks.getPositions(), blocks.m_numBlocks);
            writeArrayToFile(outFileStates, blocks.getStates(), blocks.m_numBlocks);
		}
		blocks.m_t += dt;
		counter ++;
	}
	end = clock();

	// Output parameters to file
	outFileParameters << "nx " << numBlocks << "\n";
	outFileParameters << "dt " << dt << "\n";
	outFileParameters << "tStop " << tStop << "\n";
	outFileParameters << blocks;

	// Close output files
	outFilePositions.close();
	outFileParameters.close();
        outFileStates.close();

	cout << "Ran " << counter << " integration steps for "<<numBlocks<<" blocks in " << ((double)end - (double)start)/
		CLOCKS_PER_SEC << " seconds" << endl;
	return 0;

}


void writeVectorToFile(ofstream & outFile, vector<double> &vec, int numBlocks)
{
	outFile.write(reinterpret_cast<char*>(&vec[0]), numBlocks*sizeof(double));
}

void writeArrayToFile(ofstream & outFile, double * array, int numBlocks)
{
	outFile.write(reinterpret_cast<char*>(array), numBlocks*sizeof(double));
}
