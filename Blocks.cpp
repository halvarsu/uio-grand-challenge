#include "headers/Block.h"

Block::Block(int numBlocks, double dt)
{
	dt               = dt;
	t                = 0;
	vPusher			= 4e-4;
	kPusher			= 4e6;
	k				= 2.3e6; // Stiffness between blocks
	L				= 0.14; // Physical length of block chain
	d				= L/(numBlocks-1); // Distance between blocks in block chain
	M				= 0.12;
	m				= M/numBlocks;
	eng              = sqrt(0.1)*sqrt(k*m);
	time_limit       = 0.002; // Unknown value
	mu_s             = 0.4;  // Unknown value
	mu_d             = 0.17;  // Unknown value
	f_N              = 1920/numBlocks;  // Unknown value
	k_0              = sqrt(39.2e9/f_N);  // Unknown value

}

double Block::springForce(double K, double D, double x1, double x2)
{
	return K*(x2-x1-D);
}

double Block::viscousForce(double v1, double v2)
{
	return eng*(v2-v1);
}


double Block::frictionForce(int i, double x, double v)
{
	double friction = 0;

	if (states[i]) {
		friction = -springForce(k_0, 0, start_positions[i], x);
		if (std::abs(friction) > mu_s * f_N) {
			states[i] = DYNAMIC;     // Change state
			timers[i] = t;  // Start timer
		}
	}
   
	// If the string is subsequently not attached
	if (!states[i]) {
		friction = -mu_d * f_N * sgn(v);
		timers[i] += dt;

		// Check the timer
		if (timers[i] > time_limit) {
			states[i] = STATIC;
			states[i] = true;
			start_positions[i] = x;
		}
	}
	return friction;

}
// I'm all alone in here! :(
