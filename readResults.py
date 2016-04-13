#!/usr/bin/python3
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mp
import seaborn as sns
import argparse


class FrictionAnalyser:

    def __init__(self, parameterFileName, dataFileName):
        self.parameterFileName = parameterFileName
        self.readParameters()
        self.readData(dataFileName)

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

    def readData(self, data_file):
        self.data = np.fromfile(data_file)
        nt = len(self.data) / self.nx
        self.data.resize(nt, self.nx)

    def pcolorplot(self, args):
        plt.pcolormesh(self.data - self.data,
                       cmap=mp.cm.get_cmap(args.colormap))
        plt.xlabel('Block index')
        plt.ylabel('Time step')
        plt.colorbar()
        plt.savefig("output/data.png")
        plt.show()

    def colorplot3d(self, args):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        shape = self.data.shape

        # Settings for quick tweaking
        colormap = mp.cm.get_cmap(args.colormap)
        rstep = shape[0] // 100
        cstep = shape[1] // 70

        # Plot the surface
        x = np.linspace(0, shape[1], shape[1])
        y = np.linspace(0, shape[0], shape[0])
        X, Y = np.meshgrid(x, y)
        Z = self.data - self.data[0]

        if args.plot_gradient:
            # Compute the gradient
            Gx, Gy = np.gradient(Z)
            # Normalize the magnitude
            G = (Gx**2+Gy**2)**0.5
            N = G/G.max()
            ax.plot_surface(X, Y, Z, shade=False,
                            facecolors=colormap(N), linewidth=0,
                            antialiased=True, rstride=rstep, cstride=cstep)
        else:
            ax.plot_surface(X, Y, Z, shade=True,
                            cmap=colormap, linewidth=0,
                            antialiased=True, rstride=rstep, cstride=cstep)
        # Add contours
        if args.contour:
            ax.contourf(X, Y, Z, zdir='y', offset=y[-1], cmap=colormap)
            ax.contourf(X, Y, Z, zdir='x', offset=0, cmap=colormap)

        plt.xlabel('Block index')
        plt.ylabel('Time step')
        plt.show()

    def plot(self, args):
        if args.plot_3d:
            self.colorplot3d(args)
        else:
            self.pcolorplot(args)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-3d", "--plot_3d", action="store_true",
                        help="Visualize the data in 3d")
    parser.add_argument("-c", "--contour", action="store_true",
                        help="Plot contours in 3d plot")
    parser.add_argument("-g", "--plot_gradient", action="store_true",
                        help="Overlay the gradient on the 3d surface")
    parser.add_argument("-cmap", "--colormap", help="Set the colormap",
                        choices=[m for m in mp.cm.cmap_d], default='viridis')
    return parser.parse_args()


if __name__ == "__main__":
    try:
        filenameParameters = "output/parameters.txt"
        filenamePositions = "output/positions.bin"
    except IOError:
        print("Could not read files")
        sys.exit(1)
    args = get_args()
    analyser = FrictionAnalyser(filenameParameters, filenamePositions)
    analyser.plot(args)
