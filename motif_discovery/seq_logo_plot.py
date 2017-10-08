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


    def entropy(self, probs, pi):
        """Computes the entropy"""
        if np.any(probs == 0):
            H = np.zeros(probs.shape[1])

            for j in range(probs.shape[1]):
                accum = 0.0
                for i in range(pi.shape[0]):
                    if probs[i][j] != 0.0:
                        # the sum in a sequence column of the probability of each letter times the logarithm in base 2
                        # of the probability of a letter divided by its background freq
                        accum += probs[i][j] * np.log2(probs[i][j])
                H[j] = -accum
            return H
        return np.sum(probs * np.log2(probs), axis=0).ravel()

    def relative_entropy(self, probs, pi):
        """Computes the relative entropy"""
        if np.any(probs == 0):
            RH = np.zeros(probs.shape[1])
            for j in range(probs.shape[1]):
                accum = 0.0
                for i in range(pi.shape[0]):
                    if probs[i][j] != 0.0 and pi[i] != 0.0:  # avoiding zeros
                        # the sum in a sequence column of the probability of each letter times the logarithm in base 2
                        # of the probability of a letter divided by its background freq
                        accum += probs[i][j] * np.log2(probs[i][j] / pi[i])
                RH[j] = accum
            return RH

        return np.sum(probs * np.log2(probs / pi),
                      axis=0).ravel()


    def compute_logo(self, probs, seqs, p0, aas, filename=None, method=None):
        if method is None:
            method = self.entropy
        entrop = method(probs, p0)
        contrib = probs * entrop  # contributions
        self.draw_logo(contrib, aas, filename=filename)

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