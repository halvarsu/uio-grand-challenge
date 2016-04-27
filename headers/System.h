#ifndef SYSTEM_H
#define SYSTEM_H

#define STATIC 1.0
#define DYNAMIC 0.0

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "Params.h"
#include "Block.h"

class Block;
struct connector;
class Params;
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
	double* m_positions;       // Position of each block
	double* m_velocities;      // Velocity of each block
	double* m_forces;          // Total froce on each block
	double* m_connectorForces; // Force from each connector
public:
	const int m_numBlocks  ;   // Number of blocks
    const int m_numConnectors; // Number of connectors
	double m_t       ;
    const double m_tStop   ;
    double m_dt      ;
    std::vector<Block*> m_blocks;


	friend std::ostream& operator <<(std::ostream& os, System const& system)
		{
			return os  << "nx " << system.m_numBlocks << "\n"
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
                       << "numConnectors " << system.m_numConnectors << "\n";
		}
	System(const Params & params);

	~System();
    // Trivial functions

	double* getPositions(){return m_positions;}

	double* getStates(){return m_states;}

	double* getForces(){return m_forces;}

	double* getConnectorForces(){return m_connectorForces;}

	// Non-trivial functions
    void simulate();

    void integrate();

    void copyParameters(const Params & params);

    void fillStatesArray(); // A very inefficient solution

};


#endif /* SYSTEM_H */
