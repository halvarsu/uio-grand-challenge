#ifndef BLOCK_H
#define BLOCK_H

#define STATIC 1.0
#define DYNAMIC 0.0

#include <iostream>
#include <cmath>

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
	double m_time_limit; // Crap, misaligned
	double* m_start_positions; // Even worse!
	double* m_states  ;
	double* m_timers;
	double* m_positions;
	double* m_velocities;
	double* m_forces;
public:
	const int m_numBlocks;
	double m_t       ;
	
	friend std::ostream& operator <<(std::ostream& os, Block const& blocks)
		{
			return os  << "vPusher " << blocks.m_vPusher << "\n"
					   << "kPusher " << blocks.m_kPusher << "\n"
					   << "k " << blocks.m_k << "\n"
					   << "L " << blocks.m_L << "\n"
					   << "M " << blocks.m_M << "\n"
					   << "m " << blocks.m_m << "\n"
					   << "eng " << blocks.m_eng << "\n";
		}
	Block(int numBlocks);

	~Block();

	double springForce(double K, double D, double x1, double x2);

	double viscousForce(double v1, double v2);

	double frictionForce(int i, double x, double v, double dt);

	void calculateForces(double dt);

	void integrate(double dt);

	double* getPositions();

	double* getStates();
};


#endif /* BLOCK_H */
