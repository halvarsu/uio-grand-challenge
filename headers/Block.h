#ifndef BLOCK_H
#define BLOCK_H

#include "System.h"

class System;
struct connector{
    double state;
    double x0;
    double timer;
};


class Block
{
protected:
    double m_k;                 // Spring coefficient between the blocks
    double m_m;                 // Mass of each block
    double m_f_N;               // Normal force on each block
    double m_eta;               // Dampning coefficient
    double m_connector_d;       // Distance between each connector
    double m_time_limit;        // Time for connectors to stay in dynamic
    double m_k_0;               // Spring coefficient for connectors
    double m_mu_s;              // Static friction coefficient
    double m_mu_d;              // Dynamic friction coefficient
    int m_numConnectors;        // Number of connectors/micro-junctions
    double m_dt;                // Time step
    const double* m_pT;         // Pointer to the time in system
    double m_d;                 // Distance between each block
    double m_vPusher;           // Velocity of pusher
    double m_kPusher;           // Spring coefficient of pusher
public:
    double * m_pForce;          // Pointer to total force in system
    int m_i;                    // For debugging
    double * m_pFrictionForce;  // Pointer to the friction force in system
    connector* m_connectors;    // Array of connectors
    double * m_pPosition;       // Pointer to position in system
    double * m_pVelocity;       // Pointer to velocity in system

    Block(const System& system, const int i);
    virtual ~Block();

    virtual void calculateForces();
    double springForce(const double k, const double d, const double x0, const double x1);
    double viscousForce(const double eta, const double v0, const double v1);
    double connectorForce();
    double frictionForce(const int i);
};

class FirstBlock: public Block
{
public:
    FirstBlock(System& system, const int i): Block(system, i){}
    virtual void calculateForces();
};

class LastBlock: public Block
{
public:
    LastBlock(System& system, const int i): Block(system, i){}
    virtual void calculateForces();
};

#endif /* BLOCK_H */













