#ifndef CONNECTOR_H
#define CONNECTOR_H

#include "Blocks.h"

class Connector
{
private:
    double m_startPosition;
    double m_blockPosition;
    double m_state;
    double m_timer;
    double m_k;
public:
    Connector(double k, double startPositions);
    Connector();
    ~Connector();

    /*  Trivial functions  */
    double getStartPosition(){return m_startPosition;}

    double getBlockPosition(){return m_blockPosition;}

    double getState(){return m_state;}

    double getk(){return m_k;}
    
    void setStartPosition(double x){m_startPosition = x;}

    void setBlockPosition(double x){m_blockPosition = x;}

    void setk(double k){m_k = k;}

    /* Non-trivial functions */
    int setState(double state);

    double getForce(double x, double v, double time, double dt);
};


#endif /* CONNECTOR_H */
