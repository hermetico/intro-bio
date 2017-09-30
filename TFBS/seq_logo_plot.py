import sys
from matplotlib import transforms
import matplotlib.patheffects
import matplotlib.pyplot as plt
import numpy as np
import math

COLOR_SCHEME = {'A': 'green',
                'B': 'rosybrown',
                'C': 'darkblue',
                'D': 'tan',
                'E': 'gold',
                'F': 'olivedrab',
                'G': 'orange',
                'H': 'darkgreen',
                'I': 'lightseagreen',
                'J': 'deepskyblue',
                'K': 'indigo',
                'L': 'palevioletred',
                'M': 'orange',
                'N': 'salmon',
                'O': 'purple',
                'P': 'cyan',
                'Q': 'orange',
                'R': 'yellow',
                'S': 'olive',
                'T': 'firebrick',
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


class Logo_Ploter(object):

    def create_scores(self, probs , aas):
        """Creates the scores for the plot, and sorts them"""
        scores = []
        for i in range(probs.shape[1]):
            position = []
            for j in range(aas.shape[0]):
                prob = probs[j][i]
                if probs[j][i] > 0.0:
                    position.append((aas[j].upper(),prob))
            scores.append(sorted(position, key=lambda x:x[1]))

        return scores

    def create_scores_norm(self, probs , aas):
        """Creates the scores for the plot but normalized, and sorts them"""
        scores = []
        for i in range(probs.shape[1]):
            position = []
            mmax = np.sum(probs[:, i])
            for j in range(aas.shape[0]):
                prob = probs[j][i]
                if probs[j][i] > 0.0:
                    position.append((aas[j].upper(), prob / mmax))
            scores.append(sorted(position, key=lambda x:x[1]))

        return scores

    def draw_logo_norm(self, probs, aas, size=40, filename=None):

        scores = self.create_scores_norm(probs, aas)
        num_letters = aas.shape[0]
        ymax = int(math.ceil(math.log(num_letters, 2)))

        plt.rcParams['axes.labelsize'] = 14
        plt.rcParams['axes.labelweight'] = 'bold'
        if sys.platform == 'win32':
            font_family = 'Arial'
        else:
            font_family = 'sans-serif'

        fig, ax = plt.subplots()
        fig.set_size_inches(len(scores), (ymax + 0.5) * 0.75)  ###
        ax.set_xticks(range(1, len(scores) + 1))
        ax.set_yticks(range(0, ymax + 1))
        ax.set_xticklabels(range(1, len(scores) + 1))
        ax.set_yticklabels(np.arange(0, ymax + 1))
        ax.set_ylabel('bits')
        ax.xaxis.set_tick_params(labelsize=14)
        ax.yaxis.set_tick_params(labelsize=14)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(True)

        trans_offset = transforms.offset_copy(ax.transAxes,
                                              fig=fig,
                                              x=0,
                                              y=0,
                                              units='points')
        xshift = 0
        for index, position_scores in enumerate(scores):
            max_contribution = 1. * ymax / len(position_scores)
            for base, rel_score in position_scores:
                # normalizes with maximum contribution posible on a given position
                score = rel_score * max_contribution
                txt = ax.text(0.09,
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
                yshift = window_ext.height * (score * 0.6)  ###
                trans_offset = transforms.offset_copy(txt._transform, fig=fig,
                                                      y=yshift,
                                                      units='points')
            xshift += window_ext.width * 0.8  ###
            trans_offset = transforms.offset_copy(ax.transAxes, fig=fig,
                                                  x=xshift,
                                                  units='points')
        if filename is not None:
            plt.savefig(filename)
        else:
            plt.show()

    def draw_logo(self, probs, aas, size=40, filename=None):

        scores = self.create_scores(probs, aas)
        num_letters = aas.shape[0]
        ymax = int(math.ceil(math.log(num_letters, 2)))

        plt.rcParams['axes.labelsize'] = 14
        plt.rcParams['axes.labelweight'] = 'bold'
        if sys.platform == 'win32':
            font_family = 'Arial'
        else:
            font_family = 'sans-serif'

        fig, ax = plt.subplots()
        fig.set_size_inches(len(scores), (ymax + 0.5) * 0.75)  ###
        ax.set_xticks(range(1, len(scores) + 1))
        ax.set_yticks(range(0, ymax + 1))
        ax.set_xticklabels(range(1, len(scores) + 1))
        ax.set_yticklabels(np.arange(0, ymax + 1))
        ax.set_ylabel('bits')
        ax.xaxis.set_tick_params(labelsize=14)
        ax.yaxis.set_tick_params(labelsize=14)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(True)

        trans_offset = transforms.offset_copy(ax.transAxes,
                                              fig=fig,
                                              x=0,
                                              y=0,
                                              units='points')
        xshift = 0
        for index, position_scores in enumerate(scores):
            for base, rel_score in position_scores:

                score = rel_score
                txt = ax.text(0.09,
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
                yshift = window_ext.height * (score * 0.6)  ###
                trans_offset = transforms.offset_copy(txt._transform, fig=fig,
                                                      y=yshift,
                                                      units='points')
            xshift += window_ext.width * 0.8  ###
            trans_offset = transforms.offset_copy(ax.transAxes, fig=fig,
                                                  x=xshift,
                                                  units='points')
        if filename is not None:
            plt.savefig(filename)
        else:
            plt.show()