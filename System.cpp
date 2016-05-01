#include "headers/System.h"

System::System(const Params & params): m_numBlocksX(params.m_numBlocksX),
                                       m_numBlocksY(params.m_numBlocksY),
                                       m_numConnectors(params.m_numConnectors),
                                       m_tStop(params.m_tStop),
                                       m_pusherBlockPosition(params.m_pusherBlockPosition)
{
    // Copy the rest TODO: Switch to initialization list
    copyParameters(params);

    // Compute more coefficients
	m_d				   = m_L/(m_numBlocksX-1);
	m_m				   = m_M/m_numBlocksX;
	m_eta              = sqrt(0.1)*sqrt(m_k*m_m);
    m_f_N              = m_N/m_numBlocksX;
	m_k_0              = sqrt(39.2e9*m_f_N);
    m_connector_d      = m_d/m_numConnectors;

    // Allocate arrays
	m_states		   = new double[m_numBlocksX*m_numConnectors];
	m_positions		   = new Vector[m_numBlocksY*m_numBlocksX]();
	m_velocities	   = new Vector[m_numBlocksY*m_numBlocksX]();
	m_forces		   = new Vector[m_numBlocksY*m_numBlocksX]();
	m_connectorForces  = new Vector[m_numBlocksX*m_numConnectors] ();
    // Initilize the arrays
    for (int y = 0; y < m_numBlocksY; y++) {
        // Initialize the position
        for (int x = 0; x < m_numBlocksX; x++) {
            m_positions[y*m_numBlocksX+x].x = m_d*x;
            m_positions[y*m_numBlocksX+x].y = m_d*y;
        }

    }
    createBlocks();
}

System::~System(){
    // Close all open files
    m_ofStates.close();
    m_ofForces.close();
    m_ofPositions.close();
    m_ofVelocities.close();
    m_ofConnectors.close();

    // Deallocate memory
    delete [] m_positions;
    delete [] m_velocities;
    delete [] m_forces;
    delete [] m_states;
    delete [] m_connectorForces;
}

void System::createBlocks()
{
    // Construct the block structure
    m_blocks.resize(m_numBlocksY);
    for (int y = 0; y < m_numBlocksY; y++) {
        m_blocks[y].resize(m_numBlocksX);
    }
    // Make the top layer
    for (int x = 0; x < m_numBlocksX; x++){
        m_blocks[0][x] = new TopBlock(*this, 0, x);
    }

    // Make the rest, but not the bottom
    for (int y = 1; y < m_numBlocksY-1; y++){
        for (int x = 0; x < m_numBlocksX; x++)
            m_blocks[y][x] = new Block(*this, y, x);
    }
    // Make the bottom layer
    for (int x = 0; x < m_numBlocksX; x++){
        m_blocks[m_numBlocksY-1][x] = new BottomBlock(*this, m_numBlocksY-1, x);
    }

    // Make the pusher block
    m_blocks[m_pusherBlockPosition][0] = new PusherBlock(*this, m_pusherBlockPosition, 0);

    // Link the neighbours
    /*
      TL = Top Left, L = Left, BL = Bottom Left, B = Bottom, @ = Block
      .....TL
      ..../.
      ...@-L
      ...|\.
      ...B.BL

      Format of the array:
      [TL BL L B]
     */
    for (int y = 0; y < m_numBlocksY; y++){
        for (int x = 0; x < m_numBlocksX; x++){
            if(y > 0 && x < m_numBlocksX-1) // Top left
                m_blocks[y][x]->addNeighbour(*m_blocks[y-1][x+1]);
            if(y < m_numBlocksY-1 && x < m_numBlocksX-1) // Bottom left
                m_blocks[y][x]->addNeighbour(*m_blocks[y+1][x+1]);
            if(x < m_numBlocksX-1) // Left
                m_blocks[y][x]->addNeighbour(*m_blocks[y][x+1]);
            if(y < m_numBlocksY-1) // Bottom
                m_blocks[y][x]->addNeighbour(*m_blocks[y+1][x]);
        }
    }
}

void System::copyParameters(const Params &params)
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
    m_dt         = params.m_dt;
}

void System::simulate()
{

    // Reset the and recalculate the forces
    for (int y = 0; y < m_numBlocksY; y++) {
        for (int x = 0; x < m_numBlocksX; x++) {
            m_forces[y*m_numBlocksX+x] = Vector(0,0);
        }
    }
    for (int y = 0; y < m_numBlocksY; y++) {
        for (int x = 0; x < m_numBlocksX; x++) {
            m_blocks[y][x]->calculateForces();
        }
    }
    /* Velocity-vervlet
       Algorithm:
       v_(n+1/2) = v_n + f_n*delta_t / (2*m)
       r_(n+1) = r_n + v_(n+1/2)*delta_t
       v_(n+1) = v_(n+1/2) + f_(n+1)*delta_t / (2*m)
     */
    for (int x = 0; x < m_numBlocksX; x++) {
        for (int y = 0; y < m_numBlocksY; y++) {
            Vector vel_halfstep = m_velocities[y*m_numBlocksX+x] + m_forces[y*m_numBlocksX+x]*m_dt/(2*m_m);
            m_positions[y*m_numBlocksX+x] += vel_halfstep*m_dt;
            m_velocities[y*m_numBlocksX+x] = vel_halfstep + m_forces[y*m_numBlocksX+x]*m_dt/(2*m_m);
        }
    }
}


void System::fillStatesArray()
{
    for (int x = 0; x < m_numBlocksX; x++)
    {
        for (int y = 0; y < m_numConnectors; y++) {
            m_states[x*m_numConnectors+y] = m_blocks[m_numBlocksY-1][x]->m_connectors[y].state;
        }
    }
}

void System::dumpData()
{
    fillStatesArray();
    writeArrayToFile(m_ofStates, m_states, m_numBlocksX*m_numConnectors);
    writeArrayToFile(m_ofPositions, m_positions, m_numBlocksX*m_numBlocksY);
    writeArrayToFile(m_ofVelocities, m_velocities, m_numBlocksX*m_numBlocksY);
    writeArrayToFile(m_ofForces, m_forces, m_numBlocksX*m_numBlocksY);
    writeArrayToFile(m_ofConnectors, m_connectorForces, m_numBlocksX*m_numConnectors);
}

void System::writeArrayToFile(std::ofstream & outFile, double * array, int numBlocks)
{
    if(array)
        outFile.write(reinterpret_cast<char*>(array), numBlocks*sizeof(array[0]));
}

void System::writeArrayToFile(std::ofstream & outFile, Vector * array, int numBlocks)
{
    if(array)
        outFile.write(reinterpret_cast<char*>(array), numBlocks*sizeof(Vector));
}
