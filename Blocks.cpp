#include "headers/Block.h"

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
