#include <iostream>
#include <fstream>
#include <vector>

int loadBinary(std::string filename);

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        std::cerr << "Filename needed" << std::endl;
        return 0;
    }
    if(loadBinary(std::string(argv[1])))
        return 1;
    return 0;
}

int loadBinary(std::string filename)
{
    std::ifstream infile(filename, std::ios::binary);
    if (!infile)
    {
        std::cerr << filename << " could not be opened for reading."
                  << std::endl;
        return 0;
    }

    std::vector<char> buffer((
                                 std::istreambuf_iterator<char>(infile)), 
                             (std::istreambuf_iterator<char>()));
    for(auto& c: buffer)
        std::cout << int(c) << std::endl;
}
