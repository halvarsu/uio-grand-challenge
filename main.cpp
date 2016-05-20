// Build command:
// g++ main.cpp -o main
// g++ is the compiler, main.cpp is the input file, and "-o main" specifies that the executable file is called "main".
// TODO: When the connectors reconnect, they "bounce" to a upright position
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <vector>
#include <iomanip>
#include "headers/System.h"
#include "headers/Block.h"
#include "headers/Params.h"

using namespace std;

int main() // This function runs when you execute the program.
{
	// They constructor takes care of the parameters
    Params params("params.txt");
	System system(params);
	const unsigned int writeFrequency = 1000 * system.m_tStop;
    cout << "The write frequency is " << writeFrequency << endl;

	// Create output streams. These are closed upon the deletion of System


	clock_t start, end;
	start = clock();

	unsigned int counter = 0;
    system.m_t = 0;
    system.makeNormalForce();
	while (system.m_t<system.m_tStop)
	{
        // Simulate a timestep of the entire system
		system.simulate();

		// modulo operation to check whether to write output to file on this timestep
		if ( (counter%writeFrequency) == 0){
            system.dumpData();
            std::cout << "\r" << system.m_t/system.m_tStop * 100
                      << "% completed";
            //std::cout << std::string(X, '|');
            std::cout.flush();            
        }

		system.m_t += system.m_dt;
		counter ++;
	}
    end = clock();

	// Output parameters to file
    ofstream outFileParameters("output/parameters.txt");
	outFileParameters << system;
	outFileParameters.close();

	cout << "\nRan " << counter << " integration steps for "<<system.m_numBlocksY
         <<"x" << system.m_numBlocksX << " blocks with " << system.m_numConnectors << " micro-junctions each in "
         << ((double)end - (double)start)/CLOCKS_PER_SEC << " seconds" << endl;
	return 0;

}
