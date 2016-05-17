#include "headers/System.h"
//#include <typeinfo>

System::System(const Params & params): 
                                       m_vPusher(params.m_vPusher),
                                       m_kPusher(params.m_kPusher),
                                       m_k(params.m_k),
                                       m_L(params.m_L),
                                       m_M(params.m_M),
                                       m_mu_s(params.m_mu_s),
                                       m_mu_d(params.m_mu_d),
                                       m_N(params.m_N),
                                       m_time_limit(params.m_time_limit),
                                       m_pusherBlockPosition(params.m_pusherBlockPosition),
                                       m_d(m_L/(params.m_numBlocksX-1)),
                                       m_m(m_M/params.m_numBlocksX),
                                       m_eta(sqrt(0.1)*sqrt(m_k*m_m)),
                                       m_f_N(m_N/(params.m_numBlocksX*params.m_numBlocksY)),
                                       m_k_0(sqrt(39.2e9*m_f_N)),
                                       m_connector_d(m_d/params.m_numConnectors),
                                       m_numBlocksX(params.m_numBlocksX),
                                       m_numBlocksY(params.m_numBlocksY),
                                       m_numConnectors(params.m_numConnectors),
                                       m_t(0.0),
                                       m_tStop(params.m_tStop),
                                       m_dt(params.m_dt)
{
    // Compute more coefficients


    // Open the files to be written
    openFiles(params);

    // Allocate arrays
	m_positions		   = new Vector[m_numBlocksY*m_numBlocksX]();
	m_velocities	   = new Vector[m_numBlocksY*m_numBlocksX]();
	m_forces		   = new Vector[m_numBlocksY*m_numBlocksX]();
	m_connectorForces  = new double[m_numBlocksX*m_numConnectors]();
    m_states		   = new double[m_numBlocksX*m_numConnectors];
    m_pusherForce      = new Vector[1];
    // Initilize the arrays
    for (unsigned int y = 0; y < m_numBlocksY; y++) {
        // Initialize the position
        for (unsigned int x = 0; x < m_numBlocksX; x++) {
            m_positions[y*m_numBlocksX+x].x = m_d*x;
            m_positions[y*m_numBlocksX+x].y = m_d*y;
//            std::cout << x << ',' << y << '\t' << m_forces[y*m_numBlocksX+x].x
//                      << '\t' << m_forces[y*m_numBlocksX+x].y << std::endl;
        }

    }
    createGeometry(params.m_geometry);
    linkNeighbours();
}

System::~System(){
    // Close all open files
    m_ofStates.close();
    m_ofForces.close();
    m_ofPositions.close();
    m_ofVelocities.close();
    m_ofConnectors.close();
    m_ofPusherForce.close();

    // Deallocate memory
    delete [] m_positions;
    delete [] m_velocities;
    delete [] m_forces;
    delete [] m_states;
    delete [] m_connectorForces;
    delete [] m_pusherForce;
}

int System::openFiles(const Params& params)
{
    m_ofStates.open(params.m_filenameStates);
    m_ofForces.open(params.m_filenameForces);
    m_ofPositions.open(params.m_filenamePositions);
    m_ofVelocities.open(params.m_filenameVelocities);
    m_ofConnectors.open(params.m_filenameConnectors);
    m_ofPusherForce.open(params.m_filenamePusherForce);
    return 0;
}


void System::createGeometry(const std::vector<std::vector<blockType> >& geometry)
{
    /* This works under the assumption that m_numBlocksYX is 
       set equal to that of the geometry at initializaion.
    */
    // Construct the block structure
    m_blocks.resize(m_numBlocksY);
    for (unsigned int y = 0; y < m_numBlocksY; y++) {
        m_blocks[y].resize(m_numBlocksX);
    }
    for (unsigned int y = 0; y < geometry.size(); y++){
        for (unsigned int x = 0; x < geometry[y].size(); x++){
            switch (geometry[y][x]) {
            case blockType::bBlock: {
                m_blocks[y][x] = new Block(*this, y, x);
                break;
            }
            case blockType::bTopBlock:{
                m_blocks[y][x] = new TopBlock(*this, y, x);
                break;
            }
            case blockType::bBottomBlock: {
                m_blocks[y][x] = new BottomBlock(*this, y, x);
                break;
            }
            case blockType::bPusherBlock: {
                m_blocks[y][x] = new PusherBlock(*this, y, x, m_pusherForce);
                break;
            }
            case blockType::bDeleted: {
                m_blocks[y][x] = nullptr;
                break;
            }
            default:
                std::cerr << "Unknown geometry type at x, y: "
                          << x <<", " << y << std::endl;
            }
            
        }
    }
    /*
      for (auto row: m_blocks){
      for (auto b: row){
      std::cout << typeid(b).name();
      }
      std::cout << std::endl;
      }*/
}

void System::linkNeighbours()
{
    // Link the neighbours
    /*
      TR = Top Right, R = Right, BR = Bottom Right, T = Top, @ = Block
      ...T.TR
      ...|/.
      ...@-R
      ....\.
      .....BR

      Format of the array:
      [TR BR R T]
    */
    for (unsigned int y = 0; y < m_numBlocksY; y++){
        for (unsigned int x = 0; x < m_numBlocksX; x++){
            if(m_blocks[y][x]){
                if(y > 0 && x < m_numBlocksX-1){ // Top right
                    if(m_blocks[y-1][x+1]) // Check for nullptr
                        m_blocks[y][x]->addNeighbour(*m_blocks[y-1][x+1]);
                }
                else
                    m_blocks[y][x]->setNeighbourNullptr();

                if(y < m_numBlocksY-1 && x < m_numBlocksX-1){ // Bottom right
                    if(m_blocks[y+1][x+1])
                        m_blocks[y][x]->addNeighbour(*m_blocks[y+1][x+1]);
                }
                else
                    m_blocks[y][x]->setNeighbourNullptr();

                if(x < m_numBlocksX-1){ // Right
                    if(m_blocks[y][x+1])
                        m_blocks[y][x]->addNeighbour(*m_blocks[y][x+1]);
                }
                else
                    m_blocks[y][x]->setNeighbourNullptr();

                if(y < m_numBlocksY-1){ // Top
                    if(m_blocks[y+1][x])
                        m_blocks[y][x]->addNeighbour(*m_blocks[y+1][x]);
                }
                else
                    m_blocks[y][x]->setNeighbourNullptr();
            }
        }
    }
}


void System::simulate()
{


    // Reset and recalculate the forces
    for (unsigned int y = 0; y < m_numBlocksY; y++) {
        for (unsigned int x = 0; x < m_numBlocksX; x++) {
            m_forces[y*m_numBlocksX+x] = Vector(0,0);
        }
    }


    for (unsigned int y = 0; y < m_numBlocksY; y++) {
        for (unsigned int x = 0; x < m_numBlocksX; x++) {
            if(m_blocks[y][x]) // Deleted blocks are nullptr
                m_blocks[y][x]->calculateForces();
        }
    }
    /* Velocity-vervlet
       Algorithm:
       v_(n+1/2) = v_n + f_n*delta_t / (2*m)
       r_(n+1) = r_n + v_(n+1/2)*delta_t
       v_(n+1) = v_(n+1/2) + f_(n+1)*delta_t / (2*m)
     */
    for (unsigned int y = 0; y < m_numBlocksY; y++) {
        for (unsigned int x = 0; x < m_numBlocksX; x++) {
            Vector vel_halfstep = m_velocities[y*m_numBlocksX+x] + m_forces[y*m_numBlocksX+x]*m_dt/(2*m_m);
            m_positions[y*m_numBlocksX+x] += vel_halfstep*m_dt;
            m_velocities[y*m_numBlocksX+x] = vel_halfstep + m_forces[y*m_numBlocksX+x]*m_dt/(2*m_m);
        }
    }
}


void System::fillStatesArray()
{
    for (unsigned int x = 0; x < m_numBlocksX; x++)
    {
        for (unsigned int i = 0; i < m_numConnectors; i++) {
            if(m_blocks[0][x]) // Check for nullptr
                m_states[x*m_numConnectors+i] = m_blocks[0][x]->getStateOfConnector(i);
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
    writeArrayToFile(m_ofPusherForce, m_pusherForce, 1);
}
































void System::writeArrayToFile(std::ofstream & outFile,  double * array,
const unsigned int numBlocks)
{
    if(array)
        outFile.write(reinterpret_cast<char*>(array), numBlocks*sizeof(array[0]));
}

void System::writeArrayToFile(std::ofstream & outFile, Vector * array,
const unsigned int numBlocks)
{
    if(array)
        outFile.write(reinterpret_cast<char*>(array), numBlocks*sizeof(Vector));
}
