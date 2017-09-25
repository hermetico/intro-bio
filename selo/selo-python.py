from math import log
#regular expressions lib to extract values within single quotes
import re

from matplotlib import transforms
import matplotlib as mpl
from matplotlib.font_manager import FontProperties
import matplotlib.patheffects
import matplotlib.pyplot as plt
import numpy as np
import seaborn
from gtk._gtk import PositionType

plt.style.use('seaborn-ticks')

def elog( x ):
    if x == 0.0:
        return -np.inf
    return log(x, 2)

def different_letters( lines ):
    seen = []
    for line in lines:
        for ch in line:
            if ch not in seen:
                seen.append(ch)
    return seen

def background_freq( aas, seqs ):
    freq = np.zeros((aas.shape[0]))
    for i, ch in enumerate( aas ):
        freq[i] = np.sum( seqs == ch)
    return freq #/ (seqs.shape[0] * seqs.shape[1])

def prob_distribution( aas, seqs):
    """Computes the probability of finding an aminoacid at certain 
    position in a sequence, aligning the secuences"""
    P = np.zeros((aas.shape[0], seqs.shape[1]))
    for i in range(aas.shape[0]):
        # counts all appearancespp.pprint(test)
        # sum(axis=0) sums num of appearances by column
        appearances = (seqs == aas[i]).sum(axis=0)
        P[i] += appearances
    return P * (1.0 / seqs.shape[0])

def entropy( probs, pi, seqs):
    H = np.zeros(seqs.shape[1])
    for j in range(seqs.shape[1]):
        summatory = 0.0
        for i in range(pi.shape[0]):
            if probs[i][j] != 0.0:
                summatory += probs[i][j] * log(probs[i][j], 2)
        H[j] = -summatory
    return H

def relative_entropy( probs, pi, seqs):
    RH = np.zeros(seqs.shape[1])
    for j in range(seqs.shape[1]):
        summatory = 0.0
        for i in range(pi.shape[0]):
            if probs[i][j] != 0.0 and pi[i] != 0:
                summatory += probs[i][j] * log(probs[i][j] / pi[i], 2)
        RH[j] = -summatory
    return RH

def contributions( probs, rel_ent):
    C = np.zeros(probs.shape)
    for i in range(probs.shape[0]):
        for j in range(probs.shape[1]):
            C[i][j] = probs[i][j] * rel_ent[j]
    return C
 
 
def get_from_file(fname):
    """ reads data from file and writes values into a nested list"""
    result = []
    with open(fname) as f:
        for line in f:
            innerList =  [ch for ch in line.strip()]
            if(innerList):
                result.append(innerList)
    return result
        
def create_scores( probs , aas):
    scores = []
    for i in range(probs.shape[1]):
        position = []
        for j in range(aas.shape[0]):
            if C[j][i] != 0.0:
                position.append((aas[j],C[j][i]))
        scores.append(position)

    return scores

def sort_scores(scores_to_sort):
    """sorts the scores with regards to the second value in the nested list """
    sorted_list = []
    for i in scores_to_sort:
        i = sorted(i, key = lambda x: x[1])
        sorted_list.append(i)
    return sorted_list
    
class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

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
                
        
def draw_logo(all_scores, fontfamily='Arial', size=80):
    if fontfamily == 'xkcd':
        plt.xkcd()
    else:
        mpl.rcParams['font.family'] = fontfamily

    fig, ax = plt.subplots(figsize=(len(all_scores), 2.5))

    font = FontProperties()
    font.set_size(size)
    font.set_weight('bold')
    
    font.set_family(fontfamily)

    ax.set_xticks(range(1,len(all_scores)+1))    
    ax.set_yticks(range(0,3))
    ax.set_xticklabels(range(1,len(all_scores)+1), rotation=90)
    ax.set_yticklabels(np.arange(0,3,1))    
    seaborn.despine(ax=ax, trim=True)
    
    trans_offset = transforms.offset_copy(ax.transData, 
                                          fig=fig, 
                                          x=1, 
                                          y=0, 
                                          units='dots')
   
    for index, scores in enumerate(all_scores):
        yshift = 0
        for base, score in scores:
            txt = ax.text(index+1, 
                          0, 
                          base, 
                          transform=trans_offset,
                          fontsize=80, 
                          color=COLOR_SCHEME[base],
                          ha='center',
                          fontproperties=font,

                         )
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift = window_ext.height*score
            trans_offset = transforms.offset_copy(txt._transform, 
                                                  fig=fig,
                                                  y=yshift,
                                                  units='points')
        trans_offset = transforms.offset_copy(ax.transData, 
                                              fig=fig, 
                                              x=1, 
                                              y=0, 
                                              units='points')    
    plt.show()

selo_input = get_from_file("selo_input.txt")

SEQs = np.array( selo_input, dtype=str)
AAs = np.array(different_letters( selo_input ))
L = len( SEQs[0] )
PI = background_freq( AAs, SEQs )
P = prob_distribution(AAs, SEQs)
H = entropy(P, PI, SEQs)
RH = relative_entropy(P, PI, SEQs)
C = contributions( P, RH)

scores = create_scores(C, AAs)

sorted_scores = sort_scores(scores)

draw_logo(sorted_scores)