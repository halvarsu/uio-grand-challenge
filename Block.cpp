#include "./headers/Block.h"

Block::Block(const System& system, const int i): m_k(system.m_k),
                                                 m_m(system.m_m),
                                                 m_f_N(system.m_f_N),
                                                 m_eta(system.m_eta),
                                                 m_connector_d(system.m_connector_d),
                                                 m_time_limit(system.m_time_limit),
                                                 m_k_0(system.m_k_0),
                                                 m_mu_s(system.m_mu_s),
                                                 m_mu_d(system.m_mu_d),
                                                 m_numConnectors(system.m_numConnectors),
                                                 m_dt(system.m_dt),
                                                 m_d(system.m_d),
                                                 m_vPusher(system.m_vPusher),
                                                 m_kPusher(system.m_kPusher),
                                                 m_pPosition(&system.m_positions[i]),
                                                 m_pVelocity(&system.m_velocities[i]),
                                                 m_pForce(&system.m_forces[i]),
                                                 m_pFrictionForce(&system.m_connectorForces[i*m_numConnectors]),
                                                 m_pT(&system.m_t)

{
    m_i              = i;
    m_connectors     = new connector[m_numConnectors];

    for (int j = 0; j < m_numConnectors; j++) {
        m_connectors[j].state = STATIC;
        m_connectors[j].x0 = m_d*i + m_connector_d*j;
    }
}

Block::~Block()
{
    delete[] m_connectors;
    delete m_pPosition;
    delete m_pVelocity;
    delete m_pForce;
    delete m_pFrictionForce;
    delete m_pT;
}

double Block::springForce(const double k, const double d, const double x0, const double x1)
{
    return k*(x1 - x0 - d);
}

double Block::viscousForce(const double eta, const double v0, const double v1)
{
    return eta*(v1 - v0);
}

double Block::frictionForce(const int i)
{
	if (m_connectors[i].state) {
		*(m_pFrictionForce+i) = -springForce(m_k_0, 0, m_connectors[i].x0, *m_pPosition+m_connector_d*i);
		if (std::abs(m_pFrictionForce[i]) > m_mu_s * m_f_N) {
			m_connectors[i].state = DYNAMIC;     // Change state
			m_connectors[i].timer = 0;           // Start timer
		}
	}
	// If the string is subsequently not attached
	if (!m_connectors[i].state) {
		*(m_pFrictionForce+i) = -m_mu_d * m_f_N * sgn(*m_pVelocity);
		m_connectors[i].timer += m_dt;

		// Check the timer
		if (m_connectors[i].timer > m_time_limit) {
			m_connectors[i].state = STATIC;
			m_connectors[i].x0 = *m_pPosition+m_connector_d*i;
		}
	}
    return *(m_pFrictionForce+i);
}

double Block::connectorForce()
{
    double friction = 0;
    for (int i = 0; i < m_numConnectors; i++)
        friction += frictionForce(i);

    return friction;
}

void  Block::calculateForces()
{
    *m_pForce = springForce(m_k, m_d, *m_pPosition, *(m_pPosition+1))
        -springForce(m_k, m_d, *(m_pPosition-1), *m_pPosition)
        +viscousForce(m_eta, *m_pVelocity, *(m_pVelocity+1))
        -viscousForce(m_eta, *(m_pVelocity-1), *m_pVelocity)
        +connectorForce();
}

void FirstBlock::calculateForces()
{
    double pusherPosition = m_vPusher * (*m_pT);
    *m_pForce = springForce(m_k, m_d, *m_pPosition, *(m_pPosition+1))
        + springForce(m_kPusher, 0, *m_pPosition, pusherPosition)
        + viscousForce(m_eta, *m_pVelocity, *(m_pVelocity+1))
        + connectorForce();
}

void LastBlock::calculateForces()
{
    *m_pForce = springForce(m_k, -m_d, *m_pPosition, *(m_pPosition-1))
        - viscousForce(m_eta, *m_pVelocity, *(m_pVelocity-1))
        +connectorForce();
}
