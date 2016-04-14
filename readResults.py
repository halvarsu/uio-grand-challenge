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
            exec("self.%s=%f" %(words[0],float(words[1])))

    def readData(self, data_file):
        self.data = np.fromfile(data_file)
        nt = len(self.data) / self.nx
        self.data.resize(nt, self.nx)

    def pcolorplot(self, args):
        plt.pcolormesh(self.data - self.data[0],
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
        rstep = shape[0] // args.grid_size[0]
        cstep = shape[1] // args.grid_size[1]

        # Plot the surface
        x = np.linspace(0, shape[1], shape[1])
        y = np.linspace(0, shape[0], shape[0])
        X, Y = np.meshgrid(x, y)
        Z = self.data - self.data[0]

        if args.plot_gradient:
            # Compute the gradient
            Gx, Gy = np.gradient(Z)
            # Normalize the magnitude
            G = (Gx**2 + Gy**2)**0.5
            N = G / G.max()
            surf = ax.plot_surface(X, Y, Z, shade=False, alpha=args.alpha,
                                   facecolors=colormap(N), linewidth=0,
                                   antialiased=True, rstride=rstep, cstride=cstep)
        else:
            surf = ax.plot_surface(X, Y, Z, shade=True, alpha=args.alpha,
                                   cmap=colormap, linewidth=0,
                                   antialiased=True, rstride=rstep, cstride=cstep)

        # Make a colorbar
        m = mp.cm.ScalarMappable(cmap=colormap, norm=surf.norm)
        m.set_array(G if 'G' in locals() else Z)  # The data which is colored
        plt.colorbar(m)

        # Add contours
        if args.contour:
            ax.contourf(X, Y, Z, zdir='y', offset=y[-1], cmap=colormap)
            ax.contourf(X, Y, Z, zdir='x', offset=0, cmap=colormap)
            # This line acts wierdly
            # ax.contourf(X, Y, Z, zdir='z', cmap=colormap)

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
                        choices=mp.cm.cmap_d, default='viridis')
    parser.add_argument("-G", "--grid_size", nargs=2, default=[100, 70],
                        help="Define rstride and cstride", type=int)
    parser.add_argument("-a", "--alpha", type=float, default=1,
                        help="Define the alpha level")
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
