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
public:
	const int numBlocks;
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
	double* start_positions; // Even worse!
	double* states  ;
	double* timers;
	double* positions;
	double* velocities;
	double* forces;
	
	friend std::ostream& operator <<(std::ostream& os, Block const& blocks)
		{
			return os  << "vPusher " << blocks.vPusher << "\n"
					   << "kPusher " << blocks.kPusher << "\n"
					   << "k " << blocks.k << "\n"
					   << "L " << blocks.L << "\n"
					   << "M " << blocks.M << "\n"
					   << "m " << blocks.m << "\n"
					   << "eng " << blocks.eng << "\n";
		}
	Block(int numBlocks);

	~Block();

	double springForce(double K, double D, double x1, double x2);

	double viscousForce(double v1, double v2);

	double frictionForce(int i, double x, double v, double dt);

	void calculateForces(double dt);

	void integrate(double dt);
};


#endif /* BLOCK_H */
