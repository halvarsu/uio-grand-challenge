#ifndef PARAMS_H
#define PARAMS_H

#define F_DEFAULT 0.12131415
#define I_DEFAULT -1
#define S_DEFAULT ""
#include "System.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <iterator>
#include <vector>

enum class blockType{
    bDeleted,
        bBlock,
        bTopBlock,
        bBottomBlock,
        bPusherBlock
        };
class Params
{
public:
    double m_vPusher ;    
	double m_kPusher ;    
	double m_k       ;    
	double m_L       ;    
	double m_M       ;    
	double m_mu_s    ;
	double m_mu_d    ;
    double m_N       ;
	double m_time_limit; // Crap, misaligned
	int m_numBlocksX ;
    int m_numBlocksY ; 
    double m_tStop   ;
    double m_dt      ;
    int m_numConnectors;
    int m_pusherBlockPosition;
    std::vector< std::vector<blockType> > m_geometry;
    std::string m_filenameGeometry;
    std::string m_filenamePositions;
    std::string m_filenameVelocities;
    std::string m_filenameStates;
    std::string m_filenameForces;
    std::string m_filenameConnectors;
    std::string m_filenamePusherForce;

    Params(std::string filenameParameters);

    int checkParameters();
    int readParameters(std::string filenameParameters);
    int loadGeometry(std::string filenameGeometry);

};
#endif /* PARAMS_H */
