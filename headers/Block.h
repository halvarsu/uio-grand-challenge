#ifndef BLOCK_H
#define BLOCK_H

#include "System.h"
#include "Vector.h"

class System;
class Vector;
struct connector{
    double state = STATIC;
    Vector pos0  = Vector(0,0);
    double timer = 0;
};

class Block
{
protected:
    const double m_k;                 // Spring coefficient between the blocks
    const double m_m;                 // Mass of each block
    const double m_f_N;               // Normal force on each block
    const double m_eta;               // Dampning coefficient
    const Vector m_connector_d;       // Distance between each connector
    const double m_time_limit;        // Time for connectors to stay in dynamic
    const double m_k_0;               // Spring coefficient for connectors
    const double m_mu_s;              // Static friction coefficient
    const double m_mu_d;              // Dynamic friction coefficient
    const int m_numConnectors;        // Number of connectors/micro-junctions
    const double m_dt;                // Time step
    const double *const m_pT;         // Pointer to the time in system
    const double m_d;                 // Distance between each block
    const Vector m_vPusher;           // Velocity of pusher
    const double m_kPusher;           // Spring coefficient of pusher
    Block *m_pNeighbours[4];          // Array over pointers to neighbours
    int m_neighboursCounter=0;        // Keep track of the number of neighbours
public:
    const int m_row;                  // Row
    const int m_col;                  // Column
    Vector *const m_pFrictionForce;   // Pointer to the friction force in system
    connector* m_connectors;          // Array of connectors
    Vector *const m_pPosition;        // Pointer to position in system
    Vector *const m_pVelocity;        // Pointer to velocity in system
    Vector *const m_pForce;           // Pointer to total force in system

    Block(const System& system, const int row, const int col);
    Block(const Block &obj);
    virtual ~Block();

    virtual void calculateForces();
    Vector frictionForce();
    Vector viscousForce(const double eta, const Vector v0, const Vector v1);
    Vector springForce(const double k, const double d, const Vector p0, const Vector p1);
    void addNeighbour(Block& block);
    Vector calculateNeighbourForces();
};

class PusherBlock: public Block
{
public:
    PusherBlock(const System& system, const int row, const int col): Block(system,
    row, col){}
    virtual void calculateForces();
};

class TopBlock: public Block
{
public:
    TopBlock(const System& system, const int row, const int col): Block(system, row,
    col){}
    //virtual void calculateForces();
    // Anders sa at toppblokkene skal faa en normalkraft
};
class BottomBlock: public Block
{
public:
    BottomBlock(const System& system, const int row, const int col): Block(system, row,
    col){}
    virtual void calculateForces();
};

#endif /* BLOCK_H */













