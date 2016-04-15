// Build command:
// g++ main.cpp -o main
// g++ is the compiler, main.cpp is the input file, and "-o main" specifies that the executable file is called "main".

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <vector>

#define STATIC 1.0
#define DYNAMIC 0.0

using namespace std;

struct params{
	double dt      ;
	double vPusher ;    
	double kPusher ;    
	double k       ;    
	double L       ;    
	double d       ;    
	double M       ;    
	double m       ;
	double eng     ;
	double t       ;
	double mu_s    ;
	double mu_d    ;
	double k_0     ;
	double f_N     ;
	double time_limit; // Crap, misaligned
	vector<double> start_positions; // Even worse!
	vector<double> states  ;
	vector<double> timers;
	vector<double> spring;

	
	friend ostream& operator <<(ostream& os, params const& blocks)
		{
			return os  << "vPusher " << blocks.vPusher << "\n"
					   << "kPusher " << blocks.kPusher << "\n"
					   << "k " << blocks.k << "\n"
					   << "L " << blocks.L << "\n"
					   << "M " << blocks.M << "\n"
					   << "m " << blocks.m << "\n"
					   << "eng " << blocks.eng << "\n";
		}
};

// Forward declare functions
void calculateForces(params & blocks, double * forces, double * positions, double * velocities, int numBlocks);

void integrate(double dt, double m, double * forces, double * positions, double * velocities, int numBlocks);

double springForce(double k, double d, double x1, double x2);

double viscousForce(double eng, double v1, double v2);

void writeArrayToFile(ofstream & outFile, double * array, int numBlocks);
void writeVectorToFile(ofstream & outFile, vector<double> &vec, int numBlocks);

double frictionForce(params & blocks, int i, double x, double v);
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int main() // This function runs when you execute the program.
{
	// Choose parameters

	const int numBlocks		= 70;
	const double dt	   		= 1e-7;
	const double tStop 		= 0.01;

	params blocks;
	blocks.dt               = dt;
	blocks.t                = 0;
	blocks.vPusher			= 10*4e-4;
	blocks.kPusher			= 4e6;
	blocks.k				= 2.3e6; // Stiffness between blocks
	blocks.L				= 0.14; // Physical length of block chain
	blocks.d				= blocks.L/(numBlocks-1); // Distance between blocks in block chain
	blocks.M				= 0.12;
	blocks.m				= blocks.M/numBlocks;
	blocks.eng              = sqrt(0.1)*sqrt(blocks.k*blocks.m);
	blocks.time_limit       = 0.002; // Unknown value
	blocks.mu_s             = 0.4;  // Unknown value
	blocks.mu_d             = 0.17;  // Unknown value
	blocks.f_N              = 1920/numBlocks;  // Unknown value
	blocks.k_0              = sqrt(39.2e9/blocks.f_N);  // Unknown value


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
		blocks.spring.push_back(0);
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
                        writeVectorToFile(outFileStates, blocks.spring, numBlocks);
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

void calculateForces(params & blocks, double * forces, double * positions, double * velocities, int numBlocks)
 {
	// Reset forces
	for (int i = 0; i<numBlocks; i++)
	{
		forces[i] = 0;
	}

	// First block
	double pusherPosition = blocks.vPusher*blocks.t;
	forces[0] += springForce(blocks.k, blocks.d, positions[0], positions[1])
		+ springForce(blocks.kPusher, 0, positions[0], pusherPosition)
		+ viscousForce(blocks.eng, velocities[0], velocities[1])
		+ frictionForce(blocks, 0, positions[0], velocities[0]);
        //- viscousForce(eng,velocities[0],velocities[1]);

	// Middle blocks
	for (int i = 1; i<numBlocks-1; i++)
	{
		forces[i] += springForce(blocks.k, blocks.d, positions[i], positions[i+1])
			-springForce(blocks.k, blocks.d, positions[i-1], positions[i])
			+viscousForce(blocks.eng, velocities[i], velocities[i+1])
			-viscousForce(blocks.eng, velocities[i-1],velocities[i]);
		forces[i] += frictionForce(blocks, i, positions[i], velocities[i]);
	}

	// Last block
	forces[numBlocks-1] += springForce(blocks.k, -blocks.d, positions[numBlocks-1], positions[numBlocks-2])
		- viscousForce(blocks.eng, velocities[numBlocks-1], velocities[numBlocks-2])
	    + frictionForce(blocks, numBlocks-1, positions[numBlocks-1], velocities[numBlocks-1]);
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


double springForce(double k, double d, double x1, double x2)
{
	return k*(x2-x1-d);
}
double viscousForce(double eng, double v1, double v2)
{
        return eng*(v2-v1);
}


double frictionForce(params & blocks, int i, double x, double v)
{
	double friction = 0;
	//return friction; // Remove me to activate friction
	// If the 'string' is attached, check if it is still to be attached
	if (blocks.states[i]) {
		friction = -springForce(blocks.k_0, 0, blocks.start_positions[i], x);
		if (abs(friction) > blocks.mu_s * blocks.f_N) {
			blocks.states[i] = DYNAMIC;     // Change state
			blocks.timers[i] = blocks.t;  // Start timer
		}
	}
   
	// If the string is subsequently not attached
	if (!blocks.states[i]) {
		friction = -blocks.mu_d * blocks.f_N * sgn(v);
		blocks.timers[i] += blocks.dt;

		// Check the timer
		if (blocks.timers[i] > blocks.time_limit) {
			blocks.states[i] = STATIC;
			blocks.states[i] = true;
			blocks.start_positions[i] = x;
		}
	}
	return friction;

}
void writeVectorToFile(ofstream & outFile, vector<double> &vec, int numBlocks)
{
	outFile.write(reinterpret_cast<char*>(&vec[0]), numBlocks*sizeof(double));
}

void writeArrayToFile(ofstream & outFile, double * array, int numBlocks)
{
	outFile.write(reinterpret_cast<char*>(array), numBlocks*sizeof(double));
}
