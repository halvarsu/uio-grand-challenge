#ifndef SYSTEM_H
#define SYSTEM_H

#define STATIC 1.0
#define DYNAMIC 0.0
#define SQRT2 1.4142135623730950488016887

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "Params.h"
#include "Vector.h"
#include "Block.h"

class Block;
struct connector;
class Params;
class Vector;
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class System{
private:
    friend class Block;
	double m_vPusher ;         // Velocity of the pusher
	double m_kPusher ;         // Spring coefficient of the pusher
	double m_k       ;         // Spring coefficient of the system
	double m_L       ;         // Physical length of the block chain
	double m_d       ;         // Length of a block
	double m_M       ;         // Mass of the system
	double m_m       ;         // Mass the each block
	double m_eta     ;         // Damping coefficient
	double m_mu_s    ;         // Static friction coefficient
	double m_mu_d    ;         // Dynamic friction coefficient
	double m_k_0     ;         // Spring coefficient of connectors
	double m_f_N     ;         // Normal force on each block
    double m_N       ;         // Normal force on the system
    double m_connector_d;      // Length between each connector on a block
	double m_time_limit;       // Time for connector to by in dynamic state
	double* m_states  ;        // States to be dumped to file
	Vector* m_positions;       // Position of each block
	Vector* m_velocities;      // Velocity of each block
	Vector* m_forces;          // Total froce on each block
	Vector* m_connectorForces; // Force from each connector
    const int m_pusherBlockPosition; // Index of the block pushing
    /*
      NOTE: A multidimensional array is not light weight. A quicker method is
      allocating a large block of memory, double* array = new double[sizeX*sizeY]
     */
public:
	const int m_numBlocksX  ;  // Number of blocks in the x-direction
    const int m_numBlocksY  ;  // Number of blocks in the y-direction
    const int m_numConnectors; // Number of connectors
	double m_t       ;
    const double m_tStop   ;
    double m_dt      ;
    std::vector< std::vector<Block*> > m_blocks;          // Array for pointers
                                                          // to Block objects
    std::ofstream m_ofStates;
    std::ofstream m_ofPositions;
    std::ofstream m_ofVelocities;
    std::ofstream m_ofForces;
    std::ofstream m_ofConnectors;


	friend std::ostream& operator <<(std::ostream& os, System const& system)
		{
			return os  << "nx " << system.m_numBlocksX << "\n"
                       << "ny " << system.m_numBlocksY << "\n"
                       << "dt " << system.m_dt << "\n"
                       << "tStop " << system.m_tStop << "\n"
                       << "vPusher " << system.m_vPusher << "\n"
					   << "kPusher " << system.m_kPusher << "\n"
					   << "k " << system.m_k << "\n"
					   << "L " << system.m_L << "\n"
					   << "M " << system.m_M << "\n"
					   << "m " << system.m_m << "\n"
					   << "eta " << system.m_eta << "\n"
                       << "mu_s " << system.m_mu_s << "\n"
                       << "mu_d " << system.m_mu_d << "\n"
                       << "k_0 " << system.m_k_0 << "\n"
                       << "f_N " << system.m_f_N << "\n"
                       << "N " << system.m_N << "\n"
                       << "time_limit " << system.m_time_limit << "\n"
                       << "numConnectors " << system.m_numConnectors << "\n"
                       << "pusherPosition " << system.m_pusherBlockPosition << "\n";
		}
	System(const Params & params);

	~System();
    // Trivial functions
    void openStatesFile(std::string filename){m_ofStates.open(filename);};

    void openForcesFile(std::string filename){m_ofForces.open(filename);};

    void openPositionsFile(std::string filename){m_ofPositions.open(filename);};

    void openVelocitiesFile(std::string filename){m_ofVelocities.open(filename);};

    void openConnectorsFile(std::string filename){m_ofConnectors.open(filename);};
/*	Vector* getPositions(){return m_positions;}

	double* getStates(){return m_states;}

	Vector* getForces(){return m_forces;}*/

    Vector* getConnectorForces(){return m_connectorForces;}

	// Non-trivial functions
    void createBlocks();

    void simulate();

    void integrate();

    void copyParameters(const Params & params);

    void fillStatesArray(); // A very inefficient solution

    void writeArrayToFile(std::ofstream & outFile, double * array, int numBlocks);

    void writeArrayToFile(std::ofstream & outFile, Vector * array, int numBlocks);

    void dumpData();

};


#endif /* SYSTEM_H */
