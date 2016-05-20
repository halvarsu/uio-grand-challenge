import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import SymLogNorm


class Blocks():

    def __init__(self, data, colormap):
        self.data = data
        self.fig, self.ax = plt.subplots()
        self.colormap = colormap
        self.norm = SymLogNorm(
            linthresh=0.3, linscale=0.3, vmin=-1e-8, vmax=1e-8)
        self.plot = self.ax.pcolormesh(self.data[0, :, :, 1],  # norm=self.norm,
                                       cmap=self.colormap)
        self.fig.colorbar(self.plot, ax=self.ax, extend='both')
        self.ax.set_title('Forces')
        self.ax.set_xlabel('Row index')
        self.ax.set_ylabel('Column index')
        self.min = np.amin(self.data[:, :, :, 0])
        self.max = np.amax(self.data[:, :, :, 0])
        self.animation = FuncAnimation(self.fig, self.update,
                                       interval=1, blit=True,
                                       frames=self.data.shape[0] - 1)

    def update(self, num):
        self.plot = self.ax.pcolormesh(self.data[num, :, :, 1],  # norm=self.norm,
                                       cmap=self.colormap)
        return self.plot,

    def animate(self, save=False, filename="output/animation.gif"):
        plt.show()
        if save:
            print("Writing to %s, this make take a while" % filename)
            self.animation.save(filename, writer='imagemagick', fps=30,
                                )
