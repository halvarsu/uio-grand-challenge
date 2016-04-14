// Build command:
// g++ main.cpp -o main
// g++ is the compiler, main.cpp is the input file, and "-o main" specifies that the executable file is called "main".

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>

using namespace std;

// Forward declare functions
void calculateForces(double t, double k, double d, double eng, double kPusher, double vPusher, double * forces, double * positions, double * velocities, int numBlocks);

void integrate(double dt, double m, double * forces, double * positions, double * velocities, int numBlocks);

double springForce(double k, double d, double x1, double x2);

double viscousForce(double eng, double v1, double v2);

void writeArrayToFile(ofstream & outFile, double * array, int numBlocks);

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int main() // This function runs when you execute the program.
{
	clock_t start, end;
	start = clock();
	// Choose parameters
	const int numBlocks		= 70;
	double dt			= 1e-7;
	double tStop			= 0.01;
	double t			= 0;
        

	double vPusher			= 4e-4;
	double kPusher			= 4e6;

	double k			= 2.3e6*numBlocks/70; 		
        // Stiffness between blocks (now for numBlocks amount of blocks)
	double L			= 0.14; 			
        // Physical length of block chain
	double d			= L/(numBlocks-1); 	        
        // Distance between blocks in block chain
	double M			= 0.12;
	double m			= M/numBlocks;
        double eng                      = sqrt(0.1)*sqrt(k*m);

	int writeFrequency		= 10;

	// Create output streams
	ofstream outFilePositions("output/positions.bin");
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
		positions[i] = d*i;
		velocities[i] = 0;
		forces[i] = 0;
	}

	int counter = 0;
	while (t<tStop)
	{
		// Calculate forces
		calculateForces(t, k, d, eng, kPusher, vPusher, forces, positions, velocities, numBlocks);
		integrate(dt, m, forces, positions, velocities, numBlocks);

		// modulo operation to check whether to write output to file on this timestep
		if ( (counter%writeFrequency) == 0)
		{
			writeArrayToFile(outFilePositions, positions, numBlocks);
		}
		t += dt;
		counter ++;
	}
	end = clock();

		// Output parameters to file
	outFileParameters << "nx " << numBlocks << "\n";
	outFileParameters << "dt " << dt << "\n";
	outFileParameters << "tStop " << tStop << "\n";
	outFileParameters << "vPusher " << vPusher << "\n";
	outFileParameters << "kPusher " << kPusher << "\n";
	outFileParameters << "k " << k << "\n";
	outFileParameters << "L " << L << "\n";
	outFileParameters << "M " << M << "\n";
	outFileParameters << "m " << m << "\n";
	outFileParameters << "eng " << eng << "\n";

	// Close output files
	outFilePositions.close();
	outFileParameters.close();

	cout << "Ran " << counter << " integration steps for "<<numBlocks<<" blocks in " << ((double)end - (double)start)/
		CLOCKS_PER_SEC << " seconds" << endl;
	return 0;

}

void calculateForces(double t, double k, double d, double eng, double kPusher, double vPusher, double * forces, double * positions, double * velocities, int numBlocks)
 {
	// Reset forces
	for (int i = 0; i<numBlocks; i++)
	{
		forces[i] = 0;
	}

	// First block
	double pusherPosition = vPusher*t;
	forces[0] += springForce(k, d, positions[0], positions[1]) + 
                     springForce(kPusher, 0, positions[0], pusherPosition) + 
                     viscousForce(eng, velocities[0], velocities[1]);
        //- viscousForce(eng,velocities[0],velocities[1]);

	// Middle blocks
	for (int i = 1; i<numBlocks-1; i++)
	{
		forces[i] += springForce(k, d, positions[i], positions[i+1])
			    -springForce(k, d, positions[i-1], positions[i])   
                            +viscousForce(eng, velocities[i], velocities[i+1]) 
                            -viscousForce(eng,velocities[i-1],velocities[i]);
	}

	// Last block
	forces[numBlocks-1] += springForce(k, -d, positions[numBlocks-1], positions[numBlocks-2])
                             - viscousForce(eng, velocities[numBlocks-1], velocities[numBlocks-2]);
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

void writeArrayToFile(ofstream & outFile, double * array, int numBlocks)
{
	outFile.write(reinterpret_cast<char*>(array), numBlocks*sizeof(double));
}
