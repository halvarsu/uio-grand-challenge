import numpy as np
import argparse


class Data:
    def __init__(self, filenameParameters):
        self.parameterFileName = filenameParameters
        self.readParameters()

    def readParameters(self):
        ifile = open(self.parameterFileName, "r")
        fileContents = ifile.readlines()
        for line in fileContents:
            words = line.split()
            exec("self.%s=%g" % (words[0], float(words[1])))

    def readData(self, data_file, sizeY, sizeX=0, sizeZ=2):
        data = np.fromfile(data_file)
        if sizeX == 0 and sizeZ > 1:
            nt = len(data) // int(sizeY*sizeZ)
            data.resize((nt, sizeY, sizeZ))
        elif sizeX == 0:
            nt = len(data) // int(sizeY)
            data.resize((nt, sizeY))
        else:
            # The format is [timestep, row, column, x/y]
            nt = len(data) // (int(sizeY*sizeX*sizeZ))
            data.resize((nt, sizeY, sizeX, sizeZ))
        return data

    def convert2xyz(self, inputfile, outputfile):
        def tos(vec2):
            return "("+vec2[0]+","+vec2[1]+")"
        wFile = open(outputfile, 'w')
        data = self.readData(inputfile, self.ny, self.nx)
        self.data = self.readData(inputfile, self.ny, self.nx)
        self.x = self.data[:, :, :, 0]
        self.y = self.data[:, :, :, 1]

        self.dx = self.x - np.roll(self.x, 1, axis=1)
        self.dy = self.y - np.roll(self.y, 1, axis=1)

        self.xmax = max(self.x.flatten())
        self.ymax = max(self.y.flatten())

        self.dxmax = max(self.dx.flatten())
        self.dymax = max(self.dy.flatten())

        self.x = (data[:,:,:,0]-data[0,:,:,0]+data[0,:,:,0]*10**-5).flatten() * 10**7
        self.y = (data[:,:,:,1]-data[0,:,:,1]+data[0,:,:,1]*10**-5).flatten() * 10**7
        for i in range(self.x.shape[0]):
            if i % (self.nx*self.ny) == 0:
                wFile.write(str(self.nx*self.ny)+'\nGrandC\n')
            x = str(i % self.nx)
            y = str((i//self.ny) % self.ny)
            name = tos((y, x))
            wFile.write(name + ' ' + str(self.x[i]) + ' ' + str(self.y[i]) + ' ' + '0\n')


filenameParameters = "output/parameters.txt"
filenamePositions = "output/positions.bin"
filenameXYZ = "output/positions.xyz"

if __name__ == "__main__":
    data = Data(filenameParameters)
    data.convert2xyz(filenamePositions, filenameXYZ)
