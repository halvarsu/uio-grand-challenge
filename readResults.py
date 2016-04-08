# Read binary file
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp


class FrictionAnalyser:
    def __init__(self, parameterFileName):
        self.parameterFileName = parameterFileName
        self.readParameters()

    def readParameters(self):
        ifile = open(self.parameterFileName, "r")
        fileContents = ifile.readlines()
        for line in fileContents:
            words = line.split()
            if words[0] == "nx":
                self.nx = float(words[1])
            if words[0] == "dt":
                self.dt = float(words[1])
            # TODO: get more parameters from parameter file

    def pcolorplot(self, spaceTimeFile):
        spaceTimeData = np.fromfile(spaceTimeFile)
        nt = len(spaceTimeData)/self.nx;
        spaceTimeData.resize(nt, self.nx)
        plt.pcolormesh(spaceTimeData-spaceTimeData[0])
        plt.xlabel('Block index')
        plt.ylabel('Time step')
        plt.colorbar()
        plt.savefig("output/data.png")
        plt.show()

if __name__=="__main__":
    filenameParameters = "output/parameters.txt"
    filenamePositions = "output/positions.bin"
    analyser = FrictionAnalyser(filenameParameters)
    analyser.pcolorplot(filenamePositions)
