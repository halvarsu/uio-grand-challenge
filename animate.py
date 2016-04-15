import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import cm

class Blocks():
    def __init__(self, numBlocks,L,data, point_count,colormap): 
        self.point_count     = point_count
        self.data           = data
        self.c              = data-data[0]
        self.positions      = np.linspace(0,L,numBlocks)
        self.get_color_range(colormap, self.c)

        self.fig,self.ax    = plt.subplots()
        self.scatter        = plt.scatter(self.data[0],np.zeros_like(self.data[0]), 
                                                c=self.crange[0],s=50, 
                                                cmap = colormap)

        self.animation      = FuncAnimation(self.fig, self.update,data.shape[0], 
                                                fargs=(self.data, self.scatter), 
                                                interval =25)

    def get_color_range(self,colormap, c):
        mapper = cm.ScalarMappable(cmap = colormap)
        self.crange = mapper.to_rgba(c)


    def update(self,num,data, scat):
        array           = np.array((data[num],np.zeros_like(data[num]))).T
        
        self.scatter.set_offsets(array)
        self.scatter.set_facecolor(self.crange[num]), 
        print ("%d/%d" %(num,self.point_count))
        return self.scatter,

    def animate(self, numBlocks, save = False, filename="output/animation.gif"):
        plt.show()
        if save:
            print("Writing to %s, this make take a while" %filename)
            self.animation.save(filename, writer='imagemagick', fps=30)
