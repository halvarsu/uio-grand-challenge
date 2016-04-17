#include "headers/Block.h"

Block::Block(std::string filenameParameters)
{
    readParameters(filenameParameters);
    
	m_d				   = m_L/(m_numBlocks-1); // Distance between blocks in block chain
	m_m				   = m_M/m_numBlocks;
	m_eng              = sqrt(0.1)*sqrt(m_k*m_m);
    m_f_N              = m_N/m_numBlocks;
	m_k_0              = sqrt(39.2e9/m_f_N);  

	m_start_positions  = new double[m_numBlocks];
	m_states		   = new double[m_numBlocks];
	m_positions		   = new double[m_numBlocks];
	m_timers		   = new double[m_numBlocks] ();
	m_velocities	   = new double[m_numBlocks] ();
	m_forces		   = new double[m_numBlocks] ();
	m_connectorForces  = new double[m_numBlocks] ();

	// Initialize the containers
	for (int i = 0; i < m_numBlocks; i++) {
		m_start_positions[i] = m_d*i;
		m_positions[i] = m_d*i;
		m_states[i] = STATIC;
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

double Block::springForce(double k, double d, double x1, double x2)
{
	return k*(x2-x1-d);
}

double Block::viscousForce(double v1, double v2)
{
	return m_eng*(v2-v1);
}


double Block::frictionForce(int i, double x, double v)
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
		m_timers[i] += m_dt;

		// Check the timer
		if (m_timers[i] > m_time_limit) {
			m_states[i] = STATIC;
			m_states[i] = true;
			m_start_positions[i] = x;
		}
	}
	return friction;

}

void Block::calculateForces()
{
	// Reset forces
	for (int i = 0; i<m_numBlocks; i++)
	{
		m_forces[i] = 0;
		m_connectorForces[i] = 0;
	}

	// First block
	double pusherPosition = m_vPusher*m_t;
	m_connectorForces[0] = frictionForce(0, m_positions[0], m_velocities[0]);

	m_forces[0] += springForce(m_k, m_d, m_positions[0], m_positions[1])
		+ springForce(m_kPusher, 0, m_positions[0], pusherPosition)
		+ viscousForce(m_velocities[0], m_velocities[1])
		+ m_connectorForces[0];

	// Middle blocks
	for (int i = 1; i<m_numBlocks-1; i++)
	{
		m_connectorForces[i] = frictionForce(i, m_positions[i], m_velocities[i]);

		m_forces[i] += springForce(m_k, m_d, m_positions[i], m_positions[i+1])
			-springForce(m_k, m_d, m_positions[i-1], m_positions[i])
			+viscousForce(m_velocities[i], m_velocities[i+1])
			-viscousForce(m_velocities[i-1], m_velocities[i])
			+ m_connectorForces[i];
	}

	// Last block
	m_connectorForces[m_numBlocks-1] = frictionForce(m_numBlocks-1, m_positions[m_numBlocks-1], m_velocities[m_numBlocks-1]);
	m_forces[m_numBlocks-1] = springForce(m_k, -m_d, m_positions[m_numBlocks-1], m_positions[m_numBlocks-2])
		- viscousForce(m_velocities[m_numBlocks-1], m_velocities[m_numBlocks-2])
		+ m_connectorForces[m_numBlocks-1];
}

void Block::integrate()
{
	// Euler-Cromer
	for (int i = 0; i<m_numBlocks; i++)
	{
		m_velocities[i] += m_forces[i]/m_m*m_dt;
		m_positions[i] += m_velocities[i]*m_dt;
	}
}

void Block::readParameters(std::string filenameParameters)
{
    std::ifstream infileParameters(filenameParameters);
    // ifstream is used for reading files
    if (!infileParameters)
    {
        std::cerr << "The parameter file could not be opened for reading." << std::endl;
    }

    std::string line;
    while (getline(infileParameters, line))
    {
        // Cut the string into substrings
        std::istringstream iss(line);  // Stringstream to manipulate the stream
        std::vector<std::string> tokens; // Tokens contains the substrings
        copy(std::istream_iterator<std::string>(iss), // Wierd shit
             std::istream_iterator<std::string>(),
             back_inserter(tokens));

        if (tokens[0] == "vPusher") {
            m_vPusher = atof(tokens[1].c_str());
        } else if(tokens[0] == "kPusher"){
            m_kPusher = atof(tokens[1].c_str());
        } else if(tokens[0] == "k"){
            m_k = atof(tokens[1].c_str());
        } else if(tokens[0] == "L"){
            m_L = atof(tokens[1].c_str());
        } else if(tokens[0] == "M"){
            m_M = atof(tokens[1].c_str());
        } else if(tokens[0] == "time_limit"){
            m_time_limit = atof(tokens[1].c_str());
        } else if(tokens[0] == "mu_s"){
            m_mu_s = atof(tokens[1].c_str());
        } else if(tokens[0] == "mu_d"){
            m_mu_d = atof(tokens[1].c_str());
        } else if(tokens[0] == "N"){
            m_N = atof(tokens[1].c_str());
        } else if(tokens[0] == "numBlocks"){
            m_numBlocks = atoi(tokens[1].c_str()); 
        } else if(tokens[0] == "dt"){
            m_dt = atof(tokens[1].c_str());
        } else if(tokens[0] == "tStop"){
            m_tStop = atof(tokens[1].c_str());
        }
    }
                                      
}
