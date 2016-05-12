import numpy as np
from itertools import product
from scipy.signal import convolve2d
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.artist import Artist

plt.rcParams['keymap.save'] = ''
plt.rcParams['keymap.forward'] = ''
plt.rcParams['keymap.yscale'] = ''

FILENAME = 'geometry.txt'

class Block(object):
    count_colors = ['none', 'blue', 'green', 'red', 'darkblue',
                              'darkred', 'darkgreen', 'black', 'black']
    block_color = '#DDDDDD'
    topBlock_color = '#AAAAAA'
    pusherBlock_color = 'darkred'
    bottomBlock_color = 'darkgreen'
    edge_color = '#888888'
    deleted_color = 'white'
    states = {'block': block_color, 'topBlock': topBlock_color,
              'bottomBlock': bottomBlock_color, 'pusherBlock':
              pusherBlock_color, 'deleted': deleted_color}
    valid_marks = ['selected', 'unmarked']
    
    def __init__(self, x, y, _ec=edge_color, _fc=block_color):
        self.edge_color, self.face_color = _ec, _fc
        self.poly = RegularPolygon((x + 0.5, y + 0.5),
                                   numVertices=4,
                                   radius=0.5 * np.sqrt(2),
                                   orientation=np.pi / 4,
                                   ec=_ec,
                                   fc=_fc)
        self.state = 'block'
        self.mark = 'unmarked'

    def change_state(self):
        if self.state == 'block':
            self.state = 'topBlock'
        elif self.state == 'topBlock':
            self.state = 'bottomBlock'
        elif self.state == 'bottomBlock':
            self.state = 'pusherBlock'
        elif self.state == 'pusherBlock':
            self.state = 'block'
        self.poly.set_facecolor(self.states[self.state])

    def set_state(self, state):
        if state not in self.states:
            print("Illegal state: {}".format(state))
            return

        self.state = state
        self.poly.set_facecolor(self.states[self.state])

    def mark_as(self, mark):
        if mark in self.valid_marks:
            self.mark = mark
            if mark == 'selected':
                self.poly.set_facecolor('yellow')
            elif mark == 'unmarked':
                self.poly.set_facecolor(self.states[self.state])
        else:
            print('Illegal mark: {}'.format(mark))

    def get_char(self):
        if self.state == 'block':
            return '#'
        elif self.state == 'topBlock':
            return 'T'
        elif self.state == 'bottomBlock':
            return 'B'
        elif self.state == 'pusherBlock':
            return 'P'
        elif self.state == 'deleted':
            return '-'

class Geometry(object):
    covered_color = '#DDDDDD'
    uncovered_color = '#AAAAAA'
    edge_color = '#888888'
    count_colors = ['none', 'blue', 'green', 'red', 'darkblue',
                    'darkred', 'darkgreen', 'black', 'black']
    flag_vertices = np.array([[0.25, 0.2], [0.25, 0.8],
                              [0.75, 0.65], [0.25, 0.5]])
    modes = ['normal', 'visual', 'insert', 'delete', 'modify', 'visual_line']
    mode_params = {'level': 1, 'direction': 'horizontal', 'index':0}

    def __init__(self, height, width):
        self.width, self.height = width, height

        # Create the figure and axes
        self.fig = plt.figure(figsize=((width + 2) / 3., (height + 2) / 3.))
        self.ax = self.fig.add_axes((0.05, 0.05, 0.9, 0.9),
                                    aspect='equal', frameon=False,
                                    xlim=(-0.05, width + 0.05),
                                    ylim=(-0.05, height + 0.05))

        # Create the grid of squares
        self.blocks = np.array([[Block(i, j)
                                  for j in range(height)]
                                 for i in range(width)])
        [self.ax.add_patch(sq.poly) for sq in self.blocks.flat]

        # Create event hook for mouse clicks
        self.fig.canvas.mpl_connect('button_press_event', self._button_press)
        self.fig.canvas.mpl_connect('button_release_event', self._button_release)
        self.fig.canvas.mpl_connect('key_press_event', self._key_press)
        self.fig.canvas.mpl_connect('motion_notify_event', self._motion_notify)
        #self.fig.canvas.mpl_connect('draw_event', self._draw)
        self.mode = 'normal'
        self.last_coords = None
        self.selected_blocks = []

        # Blocks to be blitted
        self.blitbuffer = []

    def draw(self, single_block=None):
        if single_block is not None:
            self.ax.draw_artist(single_block.poly)
            self.fig.canvas.blit(self.ax.bbox)
        else:
            for block in self.blitbuffer:
                self.ax.draw_artist(block.poly)
            self.fig.canvas.blit(self.ax.bbox)
            self.blitbuffer = []

    def find_closest_block(self, x, y):
        x -= 0.5;
        y -= 0.5;
        best_x = 0
        best_y = 0
        best = 99999
        for ry in range(self.height):
            for rx in range(self.width):
                dist = np.sqrt((x-rx)**2+(y-ry)**2)
                if dist <= best:
                    best, best_x, best_y = dist, rx, ry
        return best_x, best_y

    def add_to_blit(self, block):
        'The blit buffer is used to draw a batches to screen'
        self.blitbuffer.append(block)

    def _draw(self, event):
        self.background = self.fig.canvas.copy_from_bbox(self.ax.bbox)
        self.draw()

    def _click_square(self, x, y):
        self.blocks[x, y].change_state()

    def _button_press(self, event):
        if (event.xdata is None) or (event.ydata is None) or not event.inaxes:
            return

        # left mouse button: reveal square.  If the square is already clicked
        # and the correct # of mines are marked, then clear surroundig squares

        x, y = self.find_closest_block(event.xdata, event.ydata)
        if event.button == 1:
            if self.mode == 'modify':
                self._click_square(x, y)
                self.draw(self.blocks[x, y])

        # right mouse button: mark/unmark flag
        elif (event.button == 3):
            pass

        self.last_coords = np.array([x, y])

    def _button_release(self, event):
        self.last_coords = None
    
    def _key_press(self, event):
        if event.key == 'escape':
            self.set_mode('normal')
            return
        
        if self.mode == 'visual':
            self.visual_handle_key(event.key)
        elif self.mode == 'visual_line':
            self.visual_line_handle_key(event.key)
        elif self.mode == 'normal':
            if event.key == 'd':
                self.set_mode('delete')
            elif event.key == 'i':
                self.set_mode('insert')
            elif event.key == 'v':
                self.set_mode('visual')
            elif event.key == 'm':
                self.set_mode('modify')
            elif event.key == 's':
                self.save(input('Filename(ignore error message): '))
                print('Saved')
            elif event.key == 'l':
                self.set_mode('visual_line')

    def delete_block(self, block):
        """ Delete the block and remove from selected_blocks """
        block.set_state('deleted')
        if block in self.selected_blocks:
            index = np.where(self.selected_blocks == block)
            self.selected_blocks = np.delete(self.selected_blocks, index)

    def visual_handle_key(self, key):
        if key == 'd':
            for block in self.selected_blocks:
                self.delete_block(block)
                self.add_to_blit(block)
            self.draw()
        elif key == 'y':
            pass
        elif key == 'p':
            pass
        elif key == 'c':
            for block in self.selected_blocks:
                block.change_state()
                self.add_to_blit(block)
            self.draw()

    def visual_line_handle_key(self, key):
        """Like visual mode, but with lines """
        direction = self.mode_params['direction']
        limit = self.width if direction == 'vertical' else self.height
        if key == '+' and self.mode_params['index'] < limit+1:
            self.mode_params['level'] += 1
        elif key == '-' and self.mode_params['index'] > 1:
            self.mode_params['level'] -= 1
        elif key == 'h' and self.mode_params['direction'] != 'horizontal':
            self.mode_params['direction'] = 'horizontal'
            self.mode_params['level'] = 1
            self.mode_params['index'] = 0
        elif key == 'v' and self.mode_params['direction'] != 'vertical':
            self.mode_params['direction'] = 'vertical'
            self.mode_params['level'] = 1
            self.mode_params['index'] = 0
        elif key == 'd':
            for block in self.selected_blocks:
                self.delete_block(block)
                self.add_to_blit(block)
        elif key == 'c':
            for block in self.selected_blocks:
                block.change_state()
                self.add_to_blit(block)
        elif (key == 'up' or key == 'right') and self.mode_params['index'] < limit-1:
            self.mode_params['index'] += 1
        elif (key == 'down' or key == 'left') and self.mode_params['index'] > 0:
            self.mode_params['index'] -= 1
        else:
            return
        self.unselect_blocks()
        self.visual_line_build()
        self.draw()

    def visual_line_build(self):
        """Constructs the line"""
        self.selected_blocks = []
        direction = self.mode_params['direction']
        index = self.mode_params['index']
        limit = self.height if direction == 'vertical' else self.width
        for i in range(index, self.mode_params['level'] + index):
            for j in range(limit):
                if direction == 'horizontal':
                    self.select(j, i)
                else:
                    self.select(i, j)

    def select(self, y, x):
        if self.blocks[y, x].state != 'deleted':
            self.blocks[y, x].mark_as('selected')
            self.selected_blocks.append(self.blocks[y, x])
            self.add_to_blit(self.blocks[y, x])

    def _motion_notify(self, event):
        """On mouse movement"""
        x, y = event.xdata, event.ydata
        if not event.inaxes:
            return
        
        if self.last_coords is not None and self.mode == 'visual':
            """
            Mark the rectangle outlined by the mouse and append the blocks
            underneath to the self.selected_blocks
            """
            self.unselect_blocks()
            rlx, rly = self.last_coords
            rx, ry = self.find_closest_block(x, y)
            for i in np.arange(min(rlx, rx), max(rlx, rx), 1):
                for j in np.arange(min(rly, ry), max(rly, ry)+1, 1):
                    if self.blocks[i, j].mark != 'selected' and self.blocks[i, j].state!='deleted':
                        self.select(i, j)
            self.draw()
        elif self.mode == 'delete':
            try:
                rx, ry = self.find_closest_block(x, y)
                self.delete_block(self.blocks[rx, ry])
                self.draw(self.blocks[rx, ry])
            except:
                pass
            
    def set_mode(self, mode):
        'Switch to the given mode, and do the necessary clean up'
        if mode not in self.modes:
            print('Illegal mode: {}'.format(mode))
        if mode != self.mode:
            # Do the necessary changes when exiting a mode
            if self.mode == 'visual' or self.mode == 'visual_line':
                self.unselect_blocks()
                self.draw()
            
            self.mode = mode
            if mode == 'visual_line':
                self.visual_line_build()
                self.draw()

        print(self.mode)

    def unselect_blocks(self):
        for block in self.selected_blocks:
            block.mark_as('unmarked')
            self.add_to_blit(block)
        self.selected_blocks = []

    def save(self, filename):
        'Save the geometry to a txt file'
        symbols = np.zeros_like(self.blocks)
        size = self.blocks.shape
        for x in range(size[0]):
            for y in range(size[1]):
                symbols[x, y] = self.blocks[x, y].get_char()
        symbols = np.transpose(symbols)
        np.savetxt(filename, symbols, delimiter='', fmt='%s')
        
if __name__ == '__main__':
    geo = Geometry(31, 57)
    plt.show()
