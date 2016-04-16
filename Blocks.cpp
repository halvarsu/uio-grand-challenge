#include "headers/Block.h"

Block::Block(int N): m_numBlocks(N)
{
	m_t                = 0;
	m_vPusher		   = 4e-4;
	m_kPusher		   = 4e6;
	m_k			       = 2.3e6; // Stiffness between blocks
	m_L				   = 0.14; // Physical length of block chain
	m_d				   = m_L/(N-1); // Distance between blocks in block chain
	m_M				   = 0.12;
	m_m				   = m_M/N;
	m_eng              = sqrt(0.1)*sqrt(m_k*m_m);
	m_time_limit       = 0.002; 
	m_mu_s             = 0.4;  
	m_mu_d             = 0.17;  
	m_f_N              = 1920/N;  
	m_k_0              = sqrt(39.2e6/m_f_N);  

	m_start_positions  = new double[N];
	m_timers           = new double[N];
	m_positions        = new double[N];
	m_velocities       = new double[N];
	m_forces           = new double[N];
	m_states           = new double[N];
	m_connectorForces  = new double[N];

	// Initialize the containers
	for (int i = 0; i < m_numBlocks; i++) {
		m_start_positions[i] = m_d*i;			
		m_positions[i] = m_d*i;
		m_velocities[i] = 0;
		m_forces[i] = 0;
		m_states[i] = STATIC;
		m_timers[i] = 0;
		m_connectorForces[i] = 0;
	}

}

Block::~Block()
{
	delete [] m_timers;
	delete [] m_start_positions;
	delete [] m_positions;
	delete [] m_velocities;
	delete [] m_forces;
	delete [] m_states;
	delete [] m_connectorForces;
}

double Block::springForce(double K, double D, double x1, double x2)
{
	return K*(x2-x1-D);
}

double Block::viscousForce(double v1, double v2)
{
	return m_eng*(v2-v1);
}


double Block::frictionForce(int i, double x, double v, double dt)
{
	double friction = 0;

	if (m_states[i]) {
		friction = -springForce(m_k_0, 0, m_start_positions[i], x);
		if (std::abs(friction) > m_mu_s * m_f_N) {
			m_states[i] = DYNAMIC;     // Change state
			m_timers[i] = m_t;  // Start timer
		}
	}

	// If the string is subsequently not attached
	if (!m_states[i]) {
		friction = -m_mu_d * m_f_N * sgn(v);
		m_timers[i] += dt;

		// Check the timer
		if (m_timers[i] > m_time_limit) {
			m_states[i] = STATIC;
			m_states[i] = true;
			m_start_positions[i] = x;
		}
	}
	return friction;

}

void Block::calculateForces(double dt)
{
	// Reset forces
	for (int i = 0; i<m_numBlocks; i++)
	{
		m_forces[i] = 0;
		m_connectorForces[i] = 0;
	}

	// First block
	double pusherPosition = m_vPusher*m_t;
	m_connectorForces[0] = frictionForce(0, m_positions[0], m_velocities[0], dt);

	m_forces[0] += springForce(m_k, m_d, m_positions[0], m_positions[1])
		+ springForce(m_kPusher, 0, m_positions[0], pusherPosition)
		+ viscousForce(m_velocities[0], m_velocities[1])
		+ m_connectorForces[0];
	
	// Middle blocks
	for (int i = 1; i<m_numBlocks-1; i++)
	{
		m_connectorForces[i] = +frictionForce(i, m_positions[i], m_velocities[i], dt);

		m_forces[i] += springForce(m_k, m_d, m_positions[i], m_positions[i+1])
			-springForce(m_k, m_d, m_positions[i-1], m_positions[i])
			+viscousForce(m_velocities[i], m_velocities[i+1])
			-viscousForce(m_velocities[i-1], m_velocities[i])
			+ m_connectorForces[i];
	}

	// Last block
	m_connectorForces[m_numBlocks-1] = frictionForce(m_numBlocks-1, m_positions[m_numBlocks-1], m_velocities[m_numBlocks-1], dt);
	m_forces[m_numBlocks-1] = springForce(m_k, -m_d, m_positions[m_numBlocks-1], m_positions[m_numBlocks-2])
		- viscousForce(m_velocities[m_numBlocks-1], m_velocities[m_numBlocks-2])
		+ m_connectorForces[m_numBlocks-1];
}

void Block::integrate(double dt)
{
	// Euler-Cromer
	for (int i = 0; i<m_numBlocks; i++)
	{
		m_velocities[i] += m_forces[i]/m_m*dt;
		m_positions[i] += m_velocities[i]*dt;
	}
}
