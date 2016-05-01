#!/usr/bin/python3
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mp
import seaborn as sns
import argparse


class FrictionAnalyser:

    def __init__(self, filenameParameters, filenamePositions,
                 filenameStates, filenameForces, filenameConnectorForces):
        self.parameterFileName = filenameParameters
        self.readParameters()
        self.positions = self.readData(filenamePositions, self.ny, self.nx)
        self.states = self.readData(filenameStates, self.numConnectors*self.nx)
        self.forces = self.readData(filenameForces, self.ny, self.nx)
        self.connectorForces = self.readData(filenameConnectorForces, self.numConnectors*self.nx)

    def readParameters(self):
        ifile = open(self.parameterFileName, "r")
        fileContents = ifile.readlines()
        for line in fileContents:
            words = line.split()
            exec("self.%s=%f" % (words[0], float(words[1])))

    def readData(self, data_file, sizeY, sizeX=0):
        data = np.fromfile(data_file)
        if sizeX == 0:
            nt = len(data) / int(sizeY)
            data.resize((nt, sizeY))
        else:
            nt = len(data) // (int(sizeY*sizeX))
            data.resize((nt, sizeY, sizeX))
        return data

    def plot(self, args, data):
        self.plot_full(args)

    def plot_full(self, args):
        colormap = mp.cm.get_cmap(args.colormap)

        # The driving force
        forces = self.forces - self.forces[0]
        driving_force = forces[:, self.pusherPosition, 0]
        # Correct row, last column
        plt.subplot(221)
        plt.plot(range(len(driving_force)), driving_force)
        plt.xlabel('t'); plt.ylabel('F(t)')
        plt.title("Driving force")

        # The force between the connectors and the surface
        plt.subplot(222)
        plt.pcolormesh(self.connectorForces, cmap=colormap)
        plt.xlabel('Connector index'); plt.ylabel('Time step')
        plt.colorbar()
        plt.title("Connector forces")

        # The state of the connectors
        plt.subplot(223)
        #_cmap = plt.get_cmap('Greys', 2)  # Attempt at binary colobar
        plt.pcolormesh(self.states, cmap=colormap)
        plt.xlabel('Connector index')
        plt.ylabel('Time step')
        plt.title("States of the connectors")
        cax = plt.colorbar()#ticks=[0, 1])
        #cax.set_ticklabels(['Static', 'Dynamic'])
        print(self.states[0])

        # Plot the forces
        plt.subplot(224)
        plt.pcolormesh(self.forces[:, -1, :], cmap=colormap)
        plt.xlabel('Block index')
        plt.ylabel('Time step')
        plt.colorbar()
        plt.title("Forces")

        plt.tight_layout()
        plt.show()


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-cmap", "--colormap", help="Set the colormap",
                        choices=mp.cm.cmap_d, default='viridis')
    return parser.parse_args()


if __name__ == "__main__":
    filenameParameters = "output/parameters.txt"
    filenamePositions = "output/positions.bin"
    filenameStates = "output/states.bin"
    filenameForces = "output/forces.bin"
    filenameConnectorForces = "output/connectorForces.bin"

    args = get_args()
    analyser = FrictionAnalyser(filenameParameters, filenamePositions,
                                filenameStates, filenameForces,
                                filenameConnectorForces)
    analyser.plot(args, "")
