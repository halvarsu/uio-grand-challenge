#include "headers/Block.h"

Block::Block(int N): numBlocks(N)
{
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

	start_positions  = new double[N];
	timers           = new double[N];
	positions        = new double[N];
	velocities       = new double[N];
	forces           = new double[N];
	states           = new double[N];

	// Initialize the containers
	for (int i = 0; i < numBlocks; i++) {
		start_positions[i] = d*i;			
		positions[i] = d*i;
		velocities[i] = 0;
		forces[i] = 0;
		states[i] = STATIC;
		timers[i] = 0;
	}

}

Block::~Block()
{
	delete [] timers;
	delete [] start_positions;
	delete [] positions;
	delete [] velocities;
	delete [] forces;
	delete [] states;
}

double Block::springForce(double K, double D, double x1, double x2)
{
	return K*(x2-x1-D);
}

double Block::viscousForce(double v1, double v2)
{
	return eng*(v2-v1);
}


double Block::frictionForce(int i, double x, double v, double dt)
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

void Block::calculateForces(double dt)
{
	// Reset forces
	for (int i = 0; i<numBlocks; i++)
	{
		forces[i] = 0;
	}

	// First block
	double pusherPosition = vPusher*t;
	forces[0] += springForce(k, d, positions[0], positions[1])
		+ springForce(kPusher, 0, positions[0], pusherPosition)
		+ viscousForce(velocities[0], velocities[1])
		+ frictionForce(0, positions[0], velocities[0], dt);

	// Middle blocks
	for (int i = 1; i<numBlocks-1; i++)
	{
		forces[i] += springForce(k, d, positions[i], positions[i+1])
			-springForce(k, d, positions[i-1], positions[i])
			+viscousForce(velocities[i], velocities[i+1])
			-viscousForce(velocities[i-1],velocities[i])
		    +frictionForce(i, positions[i], velocities[i], dt);
	}

	// Last block
	forces[numBlocks-1] += springForce(k, -d, positions[numBlocks-1], positions[numBlocks-2])
		- viscousForce(velocities[numBlocks-1], velocities[numBlocks-2])
	    + frictionForce(numBlocks-1, positions[numBlocks-1], velocities[numBlocks-1], dt);
}

void Block::integrate(double dt)
{
	// Euler-Cromer
	for (int i = 0; i<numBlocks; i++)
	{
		velocities[i] += forces[i]/m*dt;
		positions[i] += velocities[i]*dt;
	}
}
