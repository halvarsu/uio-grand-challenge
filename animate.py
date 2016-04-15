import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Blocks():
    def __init__(self, numBlocks,L,data,colormap): 
        self.data       = data
        self.c          = data-data[0]
        self.positions  = np.linspace(0,L,numBlocks)

        self.fig,self.ax= plt.subplots()
        self.scatter    = plt.scatter(self.data[0],np.zeros_like(self.data[0]), 
                            c=self.c[-1],s=50, cmap = colormap)

        self.animation  = FuncAnimation(self.fig, self.update,data.shape[0], 
                            fargs=(self.data, self.scatter), interval =15)

    def update(self,num,data, scat):
        array           = np.array((data[num],np.zeros_like(data[num]))).T
        
        self.scatter.set_offsets(array)
        self.scatter.set_array(self.c[num]), 
        print (num)
        return self.scatter,

    def animate(self, numBlocks):
        plt.colorbar()
        plt.show()
