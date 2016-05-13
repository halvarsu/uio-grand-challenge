#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from matplotlib.widgets import Slider
import seaborn as sns
import argparse


class DiscreteSlider(Slider):
    """A matplotlib slider widget with discrete steps."""
    def __init__(self, *args, **kwargs):
        """Identical to Slider.__init__, except for the "increment" kwarg.
        "increment" specifies the step size that the slider will be discritized
        to."""
        self.inc = kwargs.pop('increment', 0.5)
        Slider.__init__(self, *args, **kwargs)

    def set_val(self, val):
        discrete_val = int(val / self.inc) * self.inc
        # We can't just call Slider.set_val(self, discrete_val), because this 
        # will prevent the slider from updating properly (it will get stuck at
        # the first step and not "slide"). Instead, we'll keep track of the
        # the continuous value as self.val and pass in the discrete value to
        # everything else.
        xy = self.poly.xy
        xy[2] = discrete_val, 1
        xy[3] = discrete_val, 0
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % discrete_val)
        if self.drawon: 
            self.ax.figure.canvas.draw()
            self.val = val
        if not self.eventson: 
            return
        for cid, func in self.observers.items():
            func(discrete_val)


class FrictionAnalyser:

    def __init__(self, filenameParameters):
        self.parameterFileName = filenameParameters
        self.readParameters()
        #self.velocities = self.readData(filenameVelocities, self.ny, self.nx)

    def readParameters(self):
        ifile = open(self.parameterFileName, "r")
        fileContents = ifile.readlines()
        for line in fileContents:
            words = line.split()
            exec("self.%s=%g" % (words[0], float(words[1])))

    def loadStates(self, filenameStates):
        self.states = self.readData(filenameStates, self.numConnectors*self.nx,
                                    sizeZ=1)
        
    def loadConnectorForces(self, filenameConnectorForces):
        self.connectorForces = self.readData(filenameConnectorForces,
                                             self.numConnectors*self.nx, sizeZ=1)

    def loadPusherForce(self, filenamePusherForce):
        self.pusherForce = self.readData(filenamePusherForce, 1)

    def loadPositions(self, filenamePositions):
        self.positions = self.readData(filenamePositions, self.ny, self.nx)

    def loadForces(self, filenameForces):
        self.forces = self.readData(filenameForces, self.ny, self.nx)

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

    def plot(self, args, data):
        self.plot_full(args)

    def plot_full(self, args):
        colormap = mp.cm.get_cmap(args.colormap)

        # The driving force
        driving_force = self.pusherForce[:, :, 0]
        fig, axes = plt.subplots(2,2)
        # Correct row, last column
        axes[0,0].plot(range(len(driving_force)), driving_force)
        axes[0][0].set_xlabel('t');
        axes[0][0].set_ylabel('F(t)')
        axes[0][0].set_title("Driving force")

        # The force between the connectors and the surface
        connectorPlot = axes[0][1].pcolormesh(self.connectorForces, cmap=colormap)
        axes[0][1].set_xlabel('Connector index');
        axes[0][1].set_ylabel('Time step')
        fig.colorbar(connectorPlot, ax=axes[0][1])
        axes[0][1].set_title("Connector forces")

        # The state of the connectors
        _cmap = plt.get_cmap('Greys', 2)  # Attempt at binary colobar
        statesPlot = axes[1][0].pcolormesh(self.states, cmap=colormap)
        axes[1][0].set_xlabel('Connector index')
        axes[1][0].set_label('Time step')
        axes[1][0].set_title("States of the connectors")
        cax = fig.colorbar(statesPlot, ax=axes[1][0])#, ticks=[0, 1])
        #cax.set_ticklabels(['Dynamic', 'Static'])

        # Plot the forces

        def update(timestep):
            forcesPlot = axes[1][1].pcolormesh(self.forces[timestep, :, :, 0], cmap=colormap)
            axes[1][1].set_xlabel('Row index')
            axes[1][1].set_ylabel('Column index')
            fig.colorbar(forcesPlot, ax=axes[1][1])
            axes[1][1].set_title("Forces")


        #sliderax = fig.add_axes([0.2, 0.02, 0.6, 0.03],
        #                             axisbg='yellow')
        #slider = DiscreteSlider(sliderax, 'Value', 0, self.forces.shape[0]-1, valinit=1, increment=1)
        #slider.on_changed(update)
        #slider.drawon = True
        #forcemagnitude = np.sqrt(forces[:, :, :, 0]**2+forces[:, :, :, 1]**2)
        #plt.pcolormesh(self.force[:, -1, :, 0], cmap=colormap)


        plt.tight_layout()
        plt.show()

    def animate(self, colormap):
        from animateforce import Blocks
        colormap = mp.cm.get_cmap(colormap)
        B = Blocks(self.forces, colormap)
        B.animate()


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-cmap", "--colormap", help="Set the colormap",
                        choices=mp.cm.cmap_d, default='viridis')
    parser.add_argument("-a", "--animate", action="store_true", default=False)
    return parser.parse_args()


if __name__ == "__main__":
    filenameParameters = "output/parameters.txt"
    filenamePositions = "output/positions.bin"
    filenameVelocities = "output/velocities.bin"
    filenameStates = "output/states.bin"
    filenameForces = "output/forces.bin"
    filenameConnectorForces = "output/connectorForces.bin"
    filenamePusherForce = "output/pusherForce.bin"

    args = get_args()
    analyser = FrictionAnalyser(filenameParameters)
    if args.animate:
        analyser.loadForces(filenameForces)
        analyser.animate(args.colormap)
    else:
        analyser.loadStates(filenameStates)
        analyser.loadConnectorForces(filenameConnectorForces)
        analyser.loadPusherForce(filenamePusherForce)
        analyser.plot(args, "")
