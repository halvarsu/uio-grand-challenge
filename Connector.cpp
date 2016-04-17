#include "headers/Connector.h"

Connector::Connector(double k, double startPosition): m_k(k)
{
    m_startPosition = startPosition;
    m_blockPosition = startPosition;
    m_state = STATIC;
    m_timer = 0;
}

Connector::Connector()
{
    m_state = STATIC;
    m_timer = 0;
}

Connector::~Connector(){}

int Connector::setState(double state){
    if(state != STATIC && state != DYNAMIC){
        std::cerr << "Invalid state!" << std::endl;
        return 0;
    }
    m_state = state;
    return 1;
}
