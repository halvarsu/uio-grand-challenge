#include "./headers/Block.h"


/*
  Constructor for Block. Initialize all of the variables from System and
  construct the array of connecturs.
 */
Block::Block(const System& system, const unsigned int row, const unsigned int col):
    m_k(system.m_k),
    m_m(system.m_m),
    m_f_N(system.m_f_N),
    m_eta(system.m_eta),
    m_time_limit(system.m_time_limit),
    m_k_0(system.m_k_0),
    m_mu_s(system.m_mu_s),
    m_mu_d(system.m_mu_d),
    m_dt(system.m_dt),
    m_d(system.m_d),
    m_vPusher(Vector(system.m_vPusher, 0)),
    m_kPusher(system.m_kPusher),
    m_pT(&system.m_t),
    m_row(row), m_col(col),
    m_pPosition(&system.m_positions[row*system.m_numBlocksX+col]),
    m_pVelocity(&system.m_velocities[row*system.m_numBlocksX+col]),
    m_pForce(&system.m_forces[row*system.m_numBlocksX+col])
{
    for (int i = 0; i < 4; i++) {
        m_pNeighbours[i] = nullptr;
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
    m_time_limit(obj.m_time_limit),
    m_k_0(obj.m_k_0),
    m_mu_s(obj.m_mu_s),
    m_mu_d(obj.m_mu_d),
    m_dt(obj.m_dt),
    m_d(obj.m_d),
    m_vPusher(obj.m_vPusher),
    m_kPusher(obj.m_kPusher),
    m_pT(obj.m_pT),
    m_neighboursCounter(0),
    m_row(obj.m_row), m_col(obj.m_col),
    m_pPosition(obj.m_pPosition),
    m_pVelocity(obj.m_pVelocity),
    m_pForce(obj.m_pForce)
{
    std::cerr << "This should never happen" << std::endl;
}
/*
Deconstructor for Block. Deallocate all pointers and array
*/
Block::~Block()
{
    delete m_pPosition;
    delete m_pVelocity;
    delete m_pForce;
    delete m_pT;
}

Vector Block::springForce(const double k, const double d, const Vector p0, const Vector p1)
{
    Vector spring = p1 - p0;
    double distance = spring.length();
    double forceTotal = (distance > 0) ? k*(distance-d)/distance: k*(distance-d);
    return forceTotal * spring;
}

double Block::springForce(const double k, const double d, const double p0, const double p1)
{
    return k*(p1 - p0 - d);
}

Vector Block::viscousForce(const double eta, const Vector v0, const Vector v1)
{
    return eta*(v1 - v0);
}

/*
  The format of the neighbours array is [top right, bottom right, right, bottom]
 */
void Block::addNeighbour(Block &block)
{
    if(m_neighboursCounter >= 4)
        std::cerr << "Can not add any more neighbours" << std::endl;
    else{
    m_pNeighbours[m_neighboursCounter] = &block;
    ++m_neighboursCounter;
    }
}

void Block::setNeighbourNullptr()
{
    if(m_neighboursCounter >= 4)
        std::cerr << "Can not add any more neighbours" << std::endl;
    else{
        m_pNeighbours[m_neighboursCounter] = nullptr;
        ++m_neighboursCounter;
    }
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
                    + viscousForce(m_eta/2, *m_pVelocity,
                                  *(neighbour->m_pVelocity));
            } else { // Orthogonal springs
                tmp = springForce(m_k, m_d, *m_pPosition,
                                  *(neighbour->m_pPosition))
                    + viscousForce(m_eta, *m_pVelocity,
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
    *m_pForce += calculateNeighbourForces();
}

PusherBlock::PusherBlock(const System& system, const unsigned int row, const unsigned int col,
                         Vector* pusherForce):
    Block(system, row, col),
    m_pPusherForce(pusherForce)
{}

void PusherBlock::calculateForces()
{
    //  
    Vector pusherPosition(m_vPusher.x * (*m_pT), m_pPosition->y);
    Vector pusherForce = springForce(m_kPusher, 0, *m_pPosition, pusherPosition);
    *m_pForce += calculateNeighbourForces()
              + pusherForce;
    *m_pPusherForce = pusherForce;
}

void BottomBlock::calculateForces()
{
    *m_pForce += calculateNeighbourForces()
    + Vector(frictionForce(), 0);
}

/*
  Constructor for BottomBlock. Initialize all of the variables from System and
  construct the array of connectors.
*/
BottomBlock::BottomBlock(const System& system, const unsigned int row, const
  unsigned int col):
    Block(system, row, col),
    m_pFrictionForce(&system.m_connectorForces[col*system.m_numConnectors]),
    m_numConnectors(system.m_numConnectors),
    m_connector_d(system.m_connector_d),
    m_dynamicLength(m_mu_d*m_f_N/m_k_0)
{
    // Create the connectors
    m_connectors = new connector[m_numConnectors];

    for (unsigned int j = 0; j < m_numConnectors; j++) {
        m_connectors[j].pos0 = m_pPosition->x;
        // Optimisation might skip struct initialisation
        m_connectors[j].state = STATIC;
        m_connectors[j].timer = 0;
    }
}

BottomBlock::~BottomBlock()
{
    //~Block();
    delete[] m_connectors;
    delete m_pFrictionForce;
}

        


double BottomBlock::frictionForce()
{
    double friction;
    for(unsigned int i = 0; i < m_numConnectors; i++)
    {
        *(m_pFrictionForce+i) = -springForce(m_k_0, 0, m_connectors[i].pos0, m_pPosition->x);
        if (m_connectors[i].state) {
            if (std::abs(m_pFrictionForce[i]) > m_mu_s * m_f_N) {
                m_connectors[i].state = DYNAMIC;     // Change state
                m_connectors[i].timer = 0;           // yStart timer
            }
        }
        // If the string is subsequently not attached
        if (!m_connectors[i].state)
        {
            double dFriction = m_mu_d*m_f_N;
            if (*(m_pFrictionForce+i) < - dFriction)
            {
                *(m_pFrictionForce+i) = - dFriction;
                m_connectors[i].pos0 = m_pPosition->x - m_dynamicLength;
            } 
            if (*(m_pFrictionForce+i) > dFriction)
            {
                *(m_pFrictionForce+i) = dFriction;
                m_connectors[i].pos0 = m_pPosition->x + m_dynamicLength;
            } 
            m_connectors[i].timer += m_dt;

            // Check the timer
            if (m_connectors[i].timer > m_time_limit) {
                m_connectors[i].state = STATIC;
                m_connectors[i].pos0 = m_pPosition->x;
            }
        }
        friction += *(m_pFrictionForce+i);
    }
    return friction;
}
/*Vector BottomBlock::frictionForce()
{
    Vector friction;
    for(int i = 0; i < m_numConnectors; i++)
    {
        if (m_connectors[i].state) {
            *(m_pFrictionForce+i) = -springForce(m_k_0, 0, m_connectors[i].pos0, m_pPosition->x+m_connector_d*i);
            if (std::abs(m_pFrictionForce[i]) > m_mu_s * m_f_N) {
                m_connectors[i].state = DYNAMIC;     // Change state
                m_connectors[i].timer = 0;           // Start timer
            }
        }
        // If the string is subsequently not attached
        if (!m_connectors[i].state) {
            *(m_pFrictionForce+i) = -m_mu_d * m_f_N * std::copysign(1.0, m_pVelocity->x);
            m_connectors[i].timer += m_dt;

            // Check the timer
            if (m_connectors[i].timer > m_time_limit) {
                m_connectors[i].state = STATIC;
                m_connectors[i].pos0 = m_pPosition->x+m_connector_d*i;
            }
        }
        friction += *(m_pFrictionForce+i);
    }
    return friction;
}

*/
