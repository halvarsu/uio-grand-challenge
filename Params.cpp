#include "headers/Params.h"

Params::Params(std::string filenameParameters)
{
    readParameters(filenameParameters);
}

void Params::readParameters(std::string filenameParameters)
{
    std::ifstream infileParameters(filenameParameters);
    // ifstream is used for reading files
    if (!infileParameters)
    {
        std::cerr << "The parameter file could not be opened for reading." << std::endl;
    }

    std::string line;
    while (getline(infileParameters, line))
    {
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
        } else if(tokens[0] == "numBlocks"){
            m_numBlocks = atoi(tokens[1].c_str()); 
        } else if(tokens[0] == "dt"){
            m_dt = atof(tokens[1].c_str());
        } else if(tokens[0] == "tStop"){
            m_tStop = atof(tokens[1].c_str());
        }
    }
    infileParameters.close();
}
