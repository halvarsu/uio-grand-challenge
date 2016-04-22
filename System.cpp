#include "headers/System.h"

System::System(Params & params): m_numBlocks(params.m_numBlocks),
    m_numConnectors(params.m_numConnectors)
{
    copyParameters(params);

	m_d				   = m_L/(m_numBlocks-1);
	m_m				   = m_M/m_numBlocks;
	m_eta              = sqrt(0.1)*sqrt(m_k*m_m);
    m_f_N              = m_N/m_numBlocks;
	m_k_0              = sqrt(39.2e9*m_f_N);
    m_connector_d      = m_d/m_numConnectors;

	m_states		   = new double[m_numBlocks*m_numConnectors];
	m_positions		   = new double[m_numBlocks];
	m_velocities	   = new double[m_numBlocks] ();
	m_forces		   = new double[m_numBlocks] ();
	m_connectorForces  = new double[m_numBlocks*m_numConnectors] ();

	// Initialize the containers
	for (int i = 0; i < m_numBlocks; i++) {
		m_positions[i] = m_d*i;
	}

    // Construct the entire block structure
    m_blocks.push_back(new FirstBlock(*this, 0));
    for (int i = 1; i < m_numBlocks-1; i++) 
        m_blocks.push_back(new Block(*this, i));
    m_blocks.push_back(new LastBlock(*this, m_numBlocks-1));
}

System::~System(){
    delete [] m_positions;
    delete [] m_velocities;
    delete [] m_forces;
    delete [] m_states;
    delete [] m_connectorForces;
}

void System::copyParameters(Params &params)
{
    m_vPusher    = params.m_vPusher;
    m_kPusher    = params.m_kPusher;
    m_k	         = params.m_k;
    m_L	         = params.m_L;
    m_d	         = params.m_d;
    m_M	         = params.m_M;
    m_mu_s	     = params.m_mu_s;
    m_mu_d	     = params.m_mu_d;
    m_k_0	     = params.m_k_0;
    m_N          = params.m_N;
    m_time_limit = params.m_time_limit;
    m_tStop      = params.m_tStop;
    m_dt         = params.m_dt;
}

void System::simulate()
{
    // Calculate the forces on each block
    // for (int i = 0; i < m_numBlocks; i++)
    //     m_blocks[i]->calculateForces();

    // // Integrate
    // // Euler-Cromer
	// for (int i = 0; i<m_numBlocks; i++)
	// {
	// 	m_velocities[i] += m_forces[i]/m_m*m_dt;
	// 	m_positions[i] += m_velocities[i]*m_dt;
	// }

    /* Velocity-vervlet
       Algorithm:
       v_(n+1/2) = v_n + f_n*delta_t / (2*m)
       r_(n+1) = r_n + v_(n+1/2)*delta_t
       v_(n+1) = v_(n+1/2) + f_(n+1)*delta_t / (2*m)
     */
    for (int i = 0; i < m_numBlocks; i++) {
        double vel_halfstep = m_velocities[i] + m_forces[i]*m_dt/(2*m_m);
        m_positions[i] += vel_halfstep*m_dt;
        m_blocks[i]->calculateForces();
        m_velocities[i] = vel_halfstep + m_forces[i]*m_dt/(2*m_m);
    }
}

void System::fillStatesArray()
{
    // Only computes for one connector
    for (int i = 0; i < m_numBlocks; i++)
    {
        for (int j = 0; j < m_numConnectors; j++) {
            m_states[i*m_numConnectors+j] = m_blocks[i]->m_connectors[j].state;
        }
    }
}
