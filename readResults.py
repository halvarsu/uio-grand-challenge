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
        self.positions = self.readData(filenamePositions, self.nx*self.ny)
        self.states = self.readData(filenameStates, self.numConnectors*self.nx)
        self.forces = self.readData(filenameForces, self.nx*self.ny)
        self.connectorForces = self.readData(filenameConnectorForces, self.numConnectors*self.nx)

    def readParameters(self):
        ifile = open(self.parameterFileName, "r")
        fileContents = ifile.readlines()
        for line in fileContents:
            words = line.split()
            exec("self.%s=%f" % (words[0], float(words[1])))

    def readData(self, data_file, resizer):
        data = np.fromfile(data_file)
        nt = len(data) // int(resizer)
        data.resize(nt, int(resizer))
        return data

    def pcolorplot(self, args, data):
        Z = data - data[0]
        print(Z.shape)
        colormap = mp.cm.get_cmap(args.colormap)
        if args.plot_gradient:
            # Compute the gradient
            Gx, Gy = np.gradient(Z)
            # Normalize the magnitude
            G = (Gx**2 + Gy**2)**0.5
            plt.pcolormesh(G, cmap=colormap)
        else:
            plt.pcolormesh(Z, cmap=colormap)
        plt.xlabel('Block index')
        plt.ylabel('Time step')
        plt.colorbar()
        plt.savefig("output/data.png")
        plt.show()

    def animate(self, args, data):
        # Animates using 50 data points
        from animate import Blocks
        colormap = mp.cm.get_cmap(args.colormap)

        point_count = 50
        data_slice = int(data.shape[0]/point_count)

        blocks = Blocks(int(self.nx), self.L, data[0::data_slice],
                        point_count, colormap)
        blocks.animate(int(self.nx), save=True)

    def colorplot3d(self, args, data):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        shape = data.shape

        # Settings for quick tweaking
        colormap = mp.cm.get_cmap(args.colormap)
        rstep = shape[0] // args.grid_size[0]
        cstep = shape[1] // args.grid_size[1]

        # Plot the surface
        x = np.linspace(0, shape[1], shape[1])
        y = np.linspace(0, shape[0], shape[0])
        X, Y = np.meshgrid(x, y)
        Z = data - data[0]

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
            ax.contourf(X, Y, Z, zdir='z', offset=0, cmap=colormap)

        plt.xlabel('Block index')
        plt.ylabel('Time step')
        plt.savefig("output/data.png")
        plt.show()

    def plot(self, args, data):
        if args.full_plot:
            self.plot_full(args)
        elif args.animate:
            self.animate(args, data)
        elif args.plot_3d:
            self.colorplot3d(args, data)
        else:
            self.pcolorplot(args, data)

    def plot_full(self, args):
        colormap = mp.cm.get_cmap(args.colormap)

        # The driving force
        forces = self.forces - self.forces[0]
        driving_force = forces[:, 0]
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
        plt.pcolormesh(self.forces, cmap=colormap)
        plt.xlabel('Block index')
        plt.ylabel('Time step')
        plt.colorbar()
        plt.title("Forces")

        plt.tight_layout()
        plt.show()


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
    parser.add_argument("-G", "--grid_size", nargs=2, default=[20, 20],
                        help="Define rstride and cstride", type=int)
    parser.add_argument("-a", "--alpha", type=float, default=1,
                        help="Define the alpha level")
    parser.add_argument("-ani", "--animate", action="store_true",
                        help="Animate the movement")
    parser.add_argument("-f", "--plot_forces", action="store_true",
                        help="Plot the forces")
    parser.add_argument("--plot_state", action="store_true",
                        help="Plots the state of each connector")
    parser.add_argument("-F", "--full_plot", action="store_true",
                        help="Plots driving force, connector forces, states and")
    parser.add_argument("--plot_connectors", action='store_true',
                        help="Plots the connector forces")
    return parser.parse_args()


if __name__ == "__main__":
    filenameParameters = "output/parameters.txt"
    filenamePositions = "output/positions.bin"
    filenameStates = "output/states.bin"
    filenameForces = "output/forces.bin"
    filenameConnectorForces = "output/connectorForces.bin"
    
    args = get_args()
    analyser = FrictionAnalyser(filenameParameters, filenamePositions,
                                filenameStates, filenameForces, filenameConnectorForces)
    if args.plot_forces:
        data = analyser.forces
    elif args.plot_state:
        data = analyser.states
    elif args.plot_connectors:
        data = analyser.connectorForces
    else:
        data = analyser.positions

    analyser.plot(args, data)
