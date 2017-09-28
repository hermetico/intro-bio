import sys
from matplotlib import transforms
import matplotlib.patheffects
import matplotlib.pyplot as plt
import numpy as np

COLOR_SCHEME = {'A': 'black',
                'B': 'rosybrown',
                'C': 'firebrick',
                'D': 'tan',
                'E': 'gold',
                'F': 'olivedrab',
                'G': 'chartreuse',
                'H': 'darkgreen',
                'I': 'lightseagreen',
                'J': 'deepskyblue',
                'K': 'blue',
                'L': 'palevioletred',
                'M': 'orange',
                'N': 'salmon',
                'O': 'purple',
                'P': 'cyan',
                'Q': 'orange',
                'R': 'yellow',
                'S': 'olive',
                'T': 'indigo',
                'U': 'dodgerblue',
                'V': 'navajowhite',
                'W': 'fuchsia',
                'X': 'brown',
                'Y': 'deeppink',
                'Z': 'darkorchid'}

class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy) + affine
        renderer.draw_path(gc, tpath, affine, rgbFace)



def draw_logo(all_scores, size=40, filename="plot.png", ymax=2):
    if sys.platform == 'win32':
        font_family = 'Arial'
    else:
        font_family = 'sans-serif'

    fig, ax = plt.subplots()
    fig.set_size_inches(len(all_scores), 5)  ###
    ax.set_xticks(range(1, len(all_scores) + 1))
    ax.set_yticks(range(0, ymax))
    ax.set_xticklabels(range(1, len(all_scores) + 1))
    ax.set_yticklabels(np.arange(0, ymax))

    trans_offset = transforms.offset_copy(ax.transAxes,
                                          fig=fig,
                                          x=0,
                                          y=0,
                                          units='points')
    xshift = 0
    for index, scores in enumerate(all_scores):
        for base, score in scores:
            txt = ax.text(0.0555,
                          0,
                          base,
                          transform=trans_offset,
                          fontsize=70,
                          color=COLOR_SCHEME[base],
                          weight='bold',
                          ha='center',
                          family=font_family)

            txt.set_clip_on(False)
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift = window_ext.height * (score * .6)  ###
            trans_offset = transforms.offset_copy(txt._transform, fig=fig,
                                                  y=yshift,
                                                  units='points')
        xshift += window_ext.width * 0.8  ###
        trans_offset = transforms.offset_copy(ax.transAxes, fig=fig,
                                              x=xshift,
                                              units='points')
    plt.savefig(filename)

    # plt.show()