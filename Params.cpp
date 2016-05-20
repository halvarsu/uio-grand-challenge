#include "headers/Params.h"


Params::Params(std::string filenameParameters):
    m_vPusher(F_DEFAULT),
    m_kPusher(F_DEFAULT),
    m_k(F_DEFAULT),
    m_L(F_DEFAULT),
    m_M(F_DEFAULT),
    m_mu_s(F_DEFAULT),
    m_mu_d(F_DEFAULT),
    m_N(F_DEFAULT),
    m_time_limit(I_DEFAULT),
    m_numBlocksX(I_DEFAULT),
    m_numBlocksY(I_DEFAULT),
    m_normalForceTime(F_DEFAULT),
    m_tStop(F_DEFAULT),
    m_dt(F_DEFAULT),
    m_numConnectors(I_DEFAULT),
    m_pusherBlockPosition(I_DEFAULT),
    m_filenameGeometry(S_DEFAULT),
    m_filenameConnectors(S_DEFAULT),
    m_filenameForces(S_DEFAULT),
    m_filenamePositions(S_DEFAULT),
    m_filenamePusherForce(S_DEFAULT),
    m_filenameStates(S_DEFAULT),
    m_filenameVelocities(S_DEFAULT)
{
    readParameters(filenameParameters);
    checkParameters();
    loadGeometry(m_filenameGeometry);
}


/*
  Read the paramters from file. It reads a text file line by line, splitting
  each line into its substrings delimited by a space. Only the first and the
  second substring are consideres, making the safe to write comments to the
  parameters. The first substring is the variable name, the second is its value.
 */
int Params::readParameters(std::string filenameParameters)
{
    std::ifstream infileParameters(filenameParameters);
    // ifstream is used for reading files
    if (!infileParameters)
    {
        std::cerr << "The parameter file could not be opened for reading." << std::endl;
        return -1;
    }

    std::string line;
    while (getline(infileParameters, line))
    {
        if(line=="") continue;
        // Cut the string into substrings
        std::istringstream iss(line);  // Stringstream to manipulate the stream
        std::vector<std::string> tokens; // Tokens contains the substrings
        copy(std::istream_iterator<std::string>(iss), // Wierd shit
             std::istream_iterator<std::string>(),
             back_inserter(tokens));

        if (tokens[0] == "vPusher") {
            m_vPusher = atof(tokens[1].c_str());
        } else if(tokens[0] == "kPusher"){
            m_kPusher = atof(tokens[1].c_str());
        } else if(tokens[0] == "k"){
            m_k = atof(tokens[1].c_str());
        } else if(tokens[0] == "L"){
            m_L = atof(tokens[1].c_str());
        } else if(tokens[0] == "M"){
            m_M = atof(tokens[1].c_str());
        } else if(tokens[0] == "time_limit"){
            m_time_limit = atof(tokens[1].c_str());
        } else if(tokens[0] == "mu_s"){
            m_mu_s = atof(tokens[1].c_str());
        } else if(tokens[0] == "mu_d"){
            m_mu_d = atof(tokens[1].c_str());
        } else if(tokens[0] == "N"){
            m_N = atof(tokens[1].c_str());
        } else if(tokens[0] == "numConnectors"){
            m_numConnectors = atoi(tokens[1].c_str()); 
        } else if(tokens[0] == "dt"){
            m_dt = atof(tokens[1].c_str());
        } else if(tokens[0] == "normalForceTime"){
            m_normalForceTime = atof(tokens[1].c_str());
        } else if(tokens[0] == "tStop"){
            m_tStop = atof(tokens[1].c_str());
        } else if(tokens[0] == "numBlocksX"){
            m_numBlocksX = atoi(tokens[1].c_str());
        } else if(tokens[0] == "numBlocksY"){
            m_numBlocksY = atoi(tokens[1].c_str());
        } else if(tokens[0] == "pusherBlockPosition"){
            m_pusherBlockPosition = atoi(tokens[1].c_str());
        } else if(tokens[0] == "filenameGeometry"){
            m_filenameGeometry = tokens[1];
        } else if(tokens[0] == "filenamePositions"){
            m_filenamePositions = tokens[1];
        } else if(tokens[0] == "filenameVelocities"){
            m_filenameVelocities = tokens[1];
        } else if(tokens[0] == "filenameStates"){
            m_filenameStates = tokens[1];
        } else if(tokens[0] == "filenameForces"){
            m_filenameForces = tokens[1];
        } else if(tokens[0] == "filenameConnectors"){
            m_filenameConnectors = tokens[1];
        } else if(tokens[0] == "filenamePusherForce"){
            m_filenamePusherForce = tokens[1];
        } 
    }
    infileParameters.close();
    return 0;
}

/*
  Make sure that all of the parameters are loaded
  If you are wondering - no, I did not type all of this out,
  but used the wonder of regex
 */
int Params::checkParameters()
{
    if(m_vPusher == F_DEFAULT){
        std::cerr << "m_vPusher is not set" << std::endl;
    }
    if(m_kPusher == F_DEFAULT){
        std::cerr << "m_kPusher is not set" << std::endl;
    }
    if(m_k == F_DEFAULT){
        std::cerr << "m_k is not set" << std::endl;
    }
    if(m_L == F_DEFAULT){
        std::cerr << "m_L is not set" << std::endl;
    }
    if(m_M == F_DEFAULT){
        std::cerr << "m_M is not set" << std::endl;
    }
    if(m_mu_s == F_DEFAULT){
        std::cerr << "m_mu_s is not set" << std::endl;
    }
    if(m_mu_d == F_DEFAULT){
        std::cerr << "m_mu_d is not set" << std::endl;
    }
    if(m_N == F_DEFAULT){
        std::cerr << "m_N is not set" << std::endl;
    }
    if(m_time_limit == I_DEFAULT){
        std::cerr << "m_time_limit is not set" << std::endl;
    }
    if(m_numBlocksX == I_DEFAULT){
        std::cerr << "m_numBlocksX is not set" << std::endl;
    }
    if(m_numBlocksY == I_DEFAULT){
        std::cerr << "m_numBlocksY is not set" << std::endl;
    }
    if(m_normalForceTime == F_DEFAULT){
        std::cerr << "m_normalForceTime is not set" << std::endl;
    }
    if(m_tStop == F_DEFAULT){
        std::cerr << "m_tStop is not set" << std::endl;
    }
    if(m_dt == F_DEFAULT){
        std::cerr << "m_dt is not set" << std::endl;
    }
    if(m_numConnectors == I_DEFAULT){
        std::cerr << "m_numConnectors is not set" << std::endl;
    }
    if(m_pusherBlockPosition == I_DEFAULT){
        std::cerr << "m_pusherBlockPosition is not set" << std::endl;
    }
    if(m_filenameGeometry == S_DEFAULT){
        std::cerr << "m_filenameGeometry is not set" << std::endl;
    }
    if(m_filenameConnectors == S_DEFAULT){
        std::cerr << "m_filenameConnectors is not set" << std::endl;
    }
    if(m_filenameForces == S_DEFAULT){
        std::cerr << "m_filenameForces is not set" << std::endl;
    }
    if(m_filenamePositions == S_DEFAULT){
        std::cerr << "m_filenamePositions is not set" << std::endl;
    }
    if(m_filenamePusherForce == S_DEFAULT){
        std::cerr << "m_filenamePusherForce is not set" << std::endl;
    }
    if(m_filenameStates == S_DEFAULT){
        std::cerr << "m_filenameStates is not set" << std::endl;
    }
    if(m_filenameVelocities == S_DEFAULT){
        std::cerr << "m_filenameVelocities is not set" << std::endl;
    }
    return 0;
}

/*
  Read a text file containing the geometry of the block structure.
  The format is an NxM "matrix" beginning at (0,0)
 */
int Params::loadGeometry(std::string filenameGeometry)
{
    std::ifstream infile(filenameGeometry);
    // ifstream is used for reading files
    if (!infile)
    {
        std::cerr << "The parameter file could not be opened for reading." << std::endl;
        return -1;
    }

    std::string line;
    while (getline(infile, line))
    {
        if(line=="") continue;
        m_geometry.push_back(std::vector<blockType>());
        for (char& c: line) {
            if (c == '#'){
                m_geometry.back().push_back(blockType::bBlock);
            }
            else if(c == 'B'){
                m_geometry.back().push_back(blockType::bBottomBlock);
            }
            else if(c == 'T'){
                m_geometry.back().push_back(blockType::bTopBlock);
            }
            else if(c == 'P'){
                m_geometry.back().push_back(blockType::bPusherBlock);
            }
            else if(c == '-'){
                m_geometry.back().push_back(blockType::bDeleted);
            }
            else
                std::cerr << "Unknown block type: " << c << std::endl;
        
        }
    }
    // check for equal length
    const unsigned int row_length = m_geometry[0].size();
    for (auto row: m_geometry){
        if (row.size() != row_length) {
            std::cerr << "Unequal length" << std::endl;
            return -1;
        }
    }
    m_numBlocksY = m_geometry.size();
    m_numBlocksX = m_geometry[0].size();
    infile.close();
    return 0;
}
