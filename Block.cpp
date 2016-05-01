#include "./headers/Block.h"


/*
  Constructor for Block. Initialize all of the variables from System and
  construct the array of connecturs.
 */
Block::Block(const System& system, const int row, const int col):
    m_k(system.m_k),
    m_m(system.m_m),
    m_f_N(system.m_f_N),
    m_eta(system.m_eta),
    m_connector_d(Vector(system.m_connector_d,0)),
    m_time_limit(system.m_time_limit),
    m_k_0(system.m_k_0),
    m_mu_s(system.m_mu_s),
    m_mu_d(system.m_mu_d),
    m_numConnectors(system.m_numConnectors),
    m_dt(system.m_dt),
    m_d(system.m_d),
    m_vPusher(Vector(system.m_vPusher, 0)),
    m_kPusher(system.m_kPusher),
    m_pPosition(&system.m_positions[row*system.m_numBlocksX+col]),
    m_pVelocity(&system.m_velocities[row*system.m_numBlocksX+col]),
    m_pForce(&system.m_forces[row*system.m_numBlocksX+col]),
    m_pFrictionForce(&system.m_connectorForces[row*m_numConnectors]),
    m_pT(&system.m_t),
    m_row(row), m_col(col)

{
    // Create the connectors
    m_connectors = new connector[m_numConnectors];
    for (int j = 0; j < m_numConnectors; j++) {
        //m_connectors[j].state = STATIC;
        m_connectors[j].pos0 = Vector(m_d*row + m_connector_d.x*j, 0);
    }
}

/*
  Copy constructor
 */
Block::Block(const Block &obj):
    m_k(obj.m_k),
    m_m(obj.m_m),
    m_f_N(obj.m_f_N),
    m_eta(obj.m_eta),
    m_connector_d(obj.m_connector_d),
    m_time_limit(obj.m_time_limit),
    m_k_0(obj.m_k_0),
    m_mu_s(obj.m_mu_s),
    m_mu_d(obj.m_mu_d),
    m_numConnectors(obj.m_numConnectors),
    m_dt(obj.m_dt),
    m_d(obj.m_d),
    m_vPusher(obj.m_vPusher),
    m_kPusher(obj.m_kPusher),
    m_pPosition(obj.m_pPosition),
    m_pVelocity(obj.m_pVelocity),
    m_pForce(obj.m_pForce),
    m_pFrictionForce(obj.m_pFrictionForce),
    m_pT(obj.m_pT),
    m_row(obj.m_row), m_col(obj.m_col),
    m_connectors(obj.m_connectors),
    m_neighboursCounter(0)
{
    std::cout<< "Yay!" << std::endl;
}
/*
Deconstructor for Block. Deallocate all pointers and array
*/
Block::~Block()
{
    delete[] m_connectors;
    delete m_pPosition;
    delete m_pVelocity;
    delete m_pForce;
    delete m_pFrictionForce;
    delete m_pT;
}

Vector Block::springForce(const double k, const double d, const Vector p0, const Vector p1)
{
    Vector spring = p1 - p0;
    double distance = spring.length();
    double forceTotal = k*(distance-d)/distance;
    return forceTotal * spring;
}

Vector Block::viscousForce(const double eta, const Vector v0, const Vector v1)
{
    return eta*(v1 - v0);
}

Vector Block::frictionForce()
{
    Vector friction;
    for(int i = 0; i < m_numConnectors; i++)
    {
        if (m_connectors[i].state) {
            *(m_pFrictionForce+i) = -springForce(m_k_0, 0, m_connectors[i].pos0, *m_pPosition+m_connector_d*i);
            if (m_pFrictionForce[i].length() > m_mu_s * m_f_N) {
                m_connectors[i].state = DYNAMIC;     // Change state
                m_connectors[i].timer = 0;           // Start timer
            }
        }
        // If the string is subsequently not attached
        if (!m_connectors[i].state) {
            *(m_pFrictionForce+i) = Vector(-m_mu_d * m_f_N * sgn(m_pVelocity->x),0);
            m_connectors[i].timer += m_dt;

            // Check the timer
            if (m_connectors[i].timer > m_time_limit) {
                m_connectors[i].state = STATIC;
                m_connectors[i].pos0 = *m_pPosition+m_connector_d*i;
            }
        }
        friction += *(m_pFrictionForce+i);
    }
    return friction;
}

/*
  The format of the neighbours array is [top left, bottom left, left, bottom]
 */
void Block::addNeighbour(Block &block)
{
    m_pNeighbours[m_neighboursCounter] = &block;
    ++m_neighboursCounter;
}

Vector Block::calculateNeighbourForces(){
    Vector force(0,0);
    for (int i = 0; i < m_neighboursCounter; i++) {
        if(m_pNeighbours[i]) // May be nullptr
        {
            Vector tmp;
            Block *neighbour = m_pNeighbours[i];
            if(i < 2){ // Diagonal springs
                tmp = springForce(m_k/2, m_d*SQRT2, *m_pPosition,
                                  *(neighbour->m_pPosition))
                    - viscousForce(m_eta/2, *m_pVelocity,
                                  *(neighbour->m_pVelocity));
            } else { // Orthogonal springs
                tmp = springForce(m_k, m_d, *m_pPosition,
                                  *(neighbour->m_pPosition))
                    - viscousForce(m_eta, *m_pVelocity,
                                  *(neighbour->m_pVelocity));
            }
            // Add the opposite force to the neighbour
            *(neighbour->m_pForce) -= tmp;

            force += tmp;
        }
    }
    return force;
}

void Block::calculateForces()
{
    *m_pForce = calculateNeighbourForces();
    /*
    *m_pForce = springForce(m_k, m_d, *m_pPosition, *(m_pPosition+1))
        -springForce(m_k, m_d, *(m_pPosition-1), *m_pPosition)
        +viscousForce(m_eta, *m_pVelocity, *(m_pVelocity+1))
        -viscousForce(m_eta, *(m_pVelocity-1), *m_pVelocity)
        +connectorForce();
    */
}

void PusherBlock::calculateForces()
{
    Vector pusherPosition = m_vPusher * (*m_pT);
    *m_pForce = calculateNeighbourForces()
              + springForce(m_kPusher, 0, *m_pPosition, pusherPosition);
    /*
    *m_pForce = springForce(m_k, m_d, *m_pPosition, *(m_pPosition+1))
        + springForce(m_kPusher, 0, *m_pPosition, pusherPosition)
        + viscousForce(m_eta, *m_pVelocity, *(m_pVelocity+1))
        + connectorForce();
    */
}

void BottomBlock::calculateForces()
{
    *m_pForce = calculateNeighbourForces()
              + frictionForce();
}
