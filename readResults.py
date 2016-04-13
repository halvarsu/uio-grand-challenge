# Read binary file
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
            elif words[0] == "dt":
                self.dt = float(words[1])
            # TODO: get more parameters from parameter file

    def pcolorplot(self, spaceTimeFile):
        spaceTimeData = np.fromfile(spaceTimeFile)
        nt = len(spaceTimeData)/self.nx;
        spaceTimeData.resize(nt, self.nx)
        plt.pcolormesh(spaceTimeData-spaceTimeData[0], cmap = "viridis")
        plt.xlabel('Block index')
        plt.ylabel('Time step')
        plt.colorbar()
        plt.savefig("output/data.png")
        plt.show()
        
    def colorplot3d(self, spaceTimeFile):
        spaceTimeData = np.fromfile(spaceTimeFile)
        nt = len(spaceTimeData)/self.nx;
        spaceTimeData.resize(nt, self.nx)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        shape =  spaceTimeData.shape 

        x = np.linspace(0,shape[1],shape[1])
        y = np.linspace(0,shape[0],shape[0])

        X,Y = np.meshgrid(x,y)
        ax.plot_surface(X,Y,spaceTimeData-spaceTimeData[0], rstride=shape[0]/20, cstride = shape[1]/20, cmap=mp.cm.viridis)
        plt.xlabel('Block index')
        plt.ylabel('Time step')
        plt.show()

if __name__=="__main__":
    filenameParameters = "output/parameters.txt"
    filenamePositions = "output/positions.bin"
    analyser = FrictionAnalyser(filenameParameters)
    try: 
        arg1 = sys.argv[1]
        if arg1 == '3d':
            analyser.colorplot3d(filenamePositions)
        else:
            print 'invalid argument, should be "3d" if any'
    except IndexError:
        analyser.pcolorplot(filenamePositions)

