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
void calculateForces(Block & blocks, double * forces, double * positions, double * velocities, int numBlocks);

void integrate(double dt, double m, double * forces, double * positions, double * velocities, int numBlocks);

void writeArrayToFile(ofstream & outFile, double * array, int numBlocks);
void writeVectorToFile(ofstream & outFile, vector<double> &vec, int numBlocks);



int main() // This function runs when you execute the program.
{
	// Choose parameters

	const int numBlocks		= 70;
	const double dt	   		= 1e-7;
	const double tStop 		= 0.01;

	// The constructor takes care of the parameters
	Block blocks(numBlocks, dt);

	int writeFrequency		= 10;

	// Create output streams
	ofstream outFilePositions("output/positions.bin");
	ofstream outFileStates("output/states.bin");
	ofstream outFileParameters("output/parameters.txt");

	// Allocate position array:
	double positions[numBlocks];
	// Allocate velocity array:
	double velocities[numBlocks];
	// Allocate force array:
	double forces[numBlocks];

	// Initialize arrays
	for (int i = 0; i<numBlocks; i++)
	{
		positions[i] = blocks.d*i;
		velocities[i] = 0;
		forces[i] = 0;
		blocks.states.push_back(STATIC); // True mean static friction
		blocks.timers.push_back(0);
		blocks.start_positions.push_back(blocks.d*i);
	}

	clock_t start, end;
	start = clock();

	int counter = 0;
	while (blocks.t<tStop)
	{
		// Calculate forces
		calculateForces(blocks, forces, positions, velocities, numBlocks);
		integrate(dt, blocks.m, forces, positions, velocities, numBlocks);

		// modulo operation to check whether to write output to file on this timestep
		if ( (counter%writeFrequency) == 0)
		{
			writeArrayToFile(outFilePositions, positions, numBlocks);
                        writeVectorToFile(outFileStates, blocks.states, numBlocks);
		}
		blocks.t += dt;
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

void calculateForces(Block & blocks, double * forces, double * positions, double * velocities, int numBlocks)
 {
	// Reset forces
	for (int i = 0; i<numBlocks; i++)
	{
		forces[i] = 0;
	}

	// First block
	double pusherPosition = blocks.vPusher*blocks.t;
	forces[0] += blocks.springForce(blocks.k, blocks.d, positions[0], positions[1])
		+ blocks.springForce(blocks.kPusher, 0, positions[0], pusherPosition)
		+ blocks.viscousForce(velocities[0], velocities[1])
		+ blocks.frictionForce(0, positions[0], velocities[0]);
        //- viscousForce(eng,velocities[0],velocities[1]);

	// Middle blocks
	for (int i = 1; i<numBlocks-1; i++)
	{
		forces[i] += blocks.springForce(blocks.k, blocks.d, positions[i], positions[i+1])
			-blocks.springForce(blocks.k, blocks.d, positions[i-1], positions[i])
			+blocks.viscousForce(velocities[i], velocities[i+1])
			-blocks.viscousForce(velocities[i-1],velocities[i]);
		forces[i] += blocks.frictionForce(i, positions[i], velocities[i]);
	}

	// Last block
	forces[numBlocks-1] += blocks.springForce(blocks.k, -blocks.d, positions[numBlocks-1],
	positions[numBlocks-2])
		- blocks.viscousForce(velocities[numBlocks-1], velocities[numBlocks-2])
	    + blocks.frictionForce(numBlocks-1, positions[numBlocks-1], velocities[numBlocks-1]);
}

void integrate(double dt, double mass, double * forces, double * positions, double * velocities, int numBlocks)
{
	// Euler-Cromer
	for (int i = 0; i<numBlocks; i++)
	{
		velocities[i] += forces[i]/mass*dt;
		positions[i] += velocities[i]*dt;
	}
}

void writeVectorToFile(ofstream & outFile, vector<double> &vec, int numBlocks)
{
	outFile.write(reinterpret_cast<char*>(&vec[0]), numBlocks*sizeof(double));
}

void writeArrayToFile(ofstream & outFile, double * array, int numBlocks)
{
	outFile.write(reinterpret_cast<char*>(array), numBlocks*sizeof(double));
}
