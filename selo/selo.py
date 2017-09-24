
#import seaborn
import matplotlib.pyplot as plt

#plt.style.use('seaborn-ticks')
from matplotlib import transforms
import matplotlib.patheffects
from matplotlib.font_manager import FontProperties
import matplotlib as mpl

import numpy as np

from math import log


def elog(x):
    if x == 0.0:
        return -np.inf
    return log(x, 2)


# In[7]:

## read file to look like this!
LINES = [
    ['G', 'A', 'H', 'O', 'D', 'E', 'F', 'S', 'E', 'X', 'R', 'C', 'C', 'K', 'S', 'C', 'S', 'G'],
    ['T', 'A', 'A', 'O', 'D', 'E', 'N', 'S', 'E', 'H', 'R', 'O', 'C', 'K', 'S', 'I', 'F', 'L'],
    ['F', 'A', 'K', 'O', 'L', 'E', 'N', 'S', 'E', 'H', 'R', 'O', 'K', 'K', 'S', 'A', 'A', 'K'],
    ['U', 'A', 'F', 'L', 'D', 'E', 'N', 'S', 'E', 'F', 'R', 'O', 'C', 'K', 'S', 'A', 'B', 'E'],
    ['D', 'H', 'A', 'O', 'D', 'E', 'E', 'S', 'E', 'K', 'R', 'O', 'C', 'K', 'S', 'A', 'C', 'E'],
    ['G', 'S', 'H', 'O', 'D', 'Y', 'N', 'S', 'E', 'K', 'R', 'O', 'C', 'K', 'S', 'Y', 'E', 'E'],
    ['O', 'H', 'H', 'O', 'D', 'E', 'N', 'S', 'E', 'L', 'R', 'O', 'C', 'K', 'S', 'T', 'E', 'J'],
    ['K', 'I', 'A', 'O', 'D', 'E', 'N', 'S', 'T', 'M', 'T', 'O', 'C', 'K', 'S', 'H', 'E', 'J']
]


def read_file(file_name):
    lines = []
    with open(file_name, 'r') as f:
        for line in f:
            # read lines and create arrays
            pass


def different_letters(lines):
    seen = []
    for line in lines:
        for ch in line:
            if ch not in seen:
                seen.append(ch)
    return seen


def background_freq(aas, seqs):
    freq = np.zeros((aas.shape[0]))
    for i, ch in enumerate(aas):
        freq[i] = np.sum(seqs == ch)
    return freq  # / (seqs.shape[0] * seqs.shape[1])


def prob_distribution(aas, seqs):
    """Computes the probability of finding an aminoacid at certain
    position in a sequence, aligning the secuences"""
    P = np.zeros((aas.shape[0], seqs.shape[1]))
    for i in range(aas.shape[0]):
        # counts all appearances
        # sum(axis=0) sums num of appearances by column
        appearances = (seqs == aas[i]).sum(axis=0)
        P[i] += appearances
    return P * (1.0 / seqs.shape[0])


def entropy(probs, pi, seqs):
    H = np.zeros(seqs.shape[1])
    for j in range(seqs.shape[1]):
        summatory = 0.0
        for i in range(pi.shape[0]):
            if probs[i][j] != 0.0:
                summatory += probs[i][j] * log(probs[i][j], 2)
        H[j] = -summatory
    return H


def relative_entropy(probs, pi, seqs):
    RH = np.zeros(seqs.shape[1])
    for j in range(seqs.shape[1]):
        summatory = 0.0
        for i in range(pi.shape[0]):
            if probs[i][j] != 0.0 and pi[i] != 0:
                summatory += probs[i][j] * log(probs[i][j] / pi[i], 2)
        RH[j] = -summatory
    return RH


def contributions(probs, rel_ent):
    C = np.zeros(probs.shape)
    for i in range(probs.shape[0]):
        for j in range(probs.shape[1]):
            C[i][j] = probs[i][j] * rel_ent[j]
    return C


def create_scores(probs, aas):
    scores = []
    for i in range(probs.shape[1]):
        position = []
        for j in range(aas.shape[0]):
            if C[j][i] != 0.0:
                position.append((aas[j], C[j][i]))
        scores.append(position)

    return scores


class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy) + affine
        renderer.draw_path(gc, tpath, affine, rgbFace)


def draw_logo(all_scores, fontfamily='Arial', size=80):
    if fontfamily == 'xkcd':
        plt.xkcd()
    else:
        mpl.rcParams['font.family'] = fontfamily

    fig, ax = plt.subplots(figsize=(len(all_scores), 2.5))

    font = FontProperties()
    font.set_size(size)
    font.set_weight('bold')

    # font.set_family(fontfamily)

    ax.set_xticks(range(1, len(all_scores) + 1))
    ax.set_yticks(range(0, 3))
    ax.set_xticklabels(range(1, len(all_scores) + 1), rotation=90)
    ax.set_yticklabels(np.arange(0, 3, 1))
    #seaborn.despine(ax=ax, trim=True)

    trans_offset = transforms.offset_copy(ax.transData,
                                          fig=fig,
                                          x=1,
                                          y=0,
                                          units='dots')

    for index, scores in enumerate(all_scores):
        yshift = 0
        for base, score in scores:
            txt = ax.text(index + 1,
                          0,
                          base,
                          transform=trans_offset,
                          fontsize=80,
                          # color=COLOR_SCHEME[base],
                          ha='center',
                          fontproperties=font,

                          )
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift = window_ext.height * score
            trans_offset = transforms.offset_copy(txt._transform,
                                                  fig=fig,
                                                  y=yshift,
                                                  units='points')
        trans_offset = transforms.offset_copy(ax.transData,
                                              fig=fig,
                                              x=1,
                                              y=0,
                                              units='points')
    #plt.show()
    fig.savefig('figure.png')


# In[8]:

## AAs or observables
SEQs = np.array(LINES, dtype=str)
AAs = np.array(different_letters(LINES))
L = len(SEQs[0])
PI = background_freq(AAs, SEQs)
P = prob_distribution(AAs, SEQs)
H = entropy(P, PI, SEQs)
RH = relative_entropy(P, PI, SEQs)
C = contributions(P, RH)

scores = create_scores(C, AAs)
draw_logo(scores)