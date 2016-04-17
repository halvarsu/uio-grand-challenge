#ifndef BLOCK_H
#define BLOCK_H

#define STATIC 1.0
#define DYNAMIC 0.0

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>
#include <stdlib.h>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class Block{
private:
	double m_vPusher ;    
	double m_kPusher ;    
	double m_k       ;    
	double m_L       ;    
	double m_d       ;    
	double m_M       ;    
	double m_m       ;
	double m_eng     ;
	double m_mu_s    ;
	double m_mu_d    ;
	double m_k_0     ;
	double m_f_N     ;
    double m_N       ;
	double m_time_limit; // Crap, misaligned
	double* m_start_positions; // Even worse!
	double* m_states  ;
	double* m_timers;
	double* m_positions;
	double* m_velocities;
	double* m_forces;
	double* m_connectorForces;
public:
	int m_numBlocks  ;
	double m_t       ;
    double m_tStop   ;
    double m_dt      ;
	
	friend std::ostream& operator <<(std::ostream& os, Block const& blocks)
		{
			return os  << "nx " << blocks.m_numBlocks << "\n"
                       << "dt " << blocks.m_dt << "\n"
                       << "tStop " << blocks.m_tStop << "\n"
                       << "vPusher " << blocks.m_vPusher << "\n"
					   << "kPusher " << blocks.m_kPusher << "\n"
					   << "k " << blocks.m_k << "\n"
					   << "L " << blocks.m_L << "\n"
					   << "M " << blocks.m_M << "\n"
					   << "m " << blocks.m_m << "\n"
					   << "eng " << blocks.m_eng << "\n"
                       << "mu_s " << blocks.m_mu_s << "\n"
                       << "mu_d " << blocks.m_mu_d << "\n"
                       << "k_0 " << blocks.m_k_0 << "\n"
                       << "f_N " << blocks.m_f_N << "\n"
                       << "N " << blocks.m_f_N << "\n"
                       << "time_limit " << blocks.m_time_limit << "\n";
		}
	Block(std::string filenameParameters);

	~Block();
    // Trivial functions

	double* getPositions(){return m_positions;}

	double* getStates(){return m_states;}

	double* getForces(){return m_forces;}

	double* getConnectorForces(){return m_connectorForces;}

	// Non-trivial functions
	double springForce(double K, double D, double x1, double x2);

	double viscousForce(double v1, double v2);

	double frictionForce(int i, double x, double v);

	void calculateForces();

	void integrate();

    void readParameters(std::string filenameParameters);
};


#endif /* BLOCK_H */
