#ifndef PARAMS_H
#define PARAMS_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <iterator>
#include <vector>


class Params
{
public:
    double m_vPusher ;    
	double m_kPusher ;    
	double m_k       ;    
	double m_L       ;    
	double m_d       ;    
	double m_M       ;    
	double m_mu_s    ;
	double m_mu_d    ;
	double m_k_0     ;
    double m_N       ;
	double m_time_limit; // Crap, misaligned
	int m_numBlocksX ;
    int m_numBlocksY ; 
    double m_tStop   ;
    double m_dt      ;
    int m_numConnectors;
    int m_pusherBlockPosition;

    Params(std::string filenameParameters);

    void readParameters(std::string filenameParameters);

};
#endif /* PARAMS_H */
