from math import log, ceil

import sys
from matplotlib import transforms
import matplotlib as mpl
from matplotlib.font_manager import FontProperties
import matplotlib.patheffects
import matplotlib.pyplot as plt
import numpy as np
#import seaborn


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

plt.style.use('seaborn-ticks')


def different_letters( lines ):
    """Returns an array with the different letters contained by the array of lines"""
    seen = {}
    for line in lines:
        for ch in line:
            if ch not in seen:
                seen[ch] = None
    return seen.keys()


def background_freq( aas, seqs ):
    """Computes the background freq"""
    freq = np.zeros((aas.shape[0]))
    for i, ch in enumerate( aas ):
        freq[i] = np.sum( seqs == ch)
    return freq


def prob_distribution( aas, seqs):
    """Computes the probability of finding an aminoacid at certain 
    position in a sequence, aligning the secuences"""
    P = np.zeros((aas.shape[0], seqs.shape[1]))
    for i in range(aas.shape[0]):
        # counts all appearances  by column
        # sum(axis=0) sums num of appearances by column
        appearances = (seqs == aas[i]).sum(axis=0)
        P[i] += appearances
    return P * (1.0 / seqs.shape[0])


def entropy( probs, pi, seqs):
    H = np.zeros(seqs.shape[1])
    for j in range(seqs.shape[1]):
        accum = 0.0
        for i in range(pi.shape[0]):
            if probs[i][j] != 0.0:
                # the sum in a sequence column of the probability of each letter times the logarithm in base 2
                # of the probability of a letter divided by its background freq
                accum += probs[i][j] * log(probs[i][j], 2)
        H[j] = -accum
    return H


def relative_entropy( probs, pi, seqs):
    RH = np.zeros(seqs.shape[1])
    for j in range(seqs.shape[1]):
        accum = 0.0
        for i in range(pi.shape[0]):
            if probs[i][j] != 0.0 and pi[i] != 0: # avoiding zeros
                # the sum in a sequence column of the probability of each letter times the logarithm in base 2
                # of the probability of a letter divided by its background freq
                accum += probs[i][j] * log(probs[i][j] / pi[i], 2)
        RH[j] = -accum
    return RH

def contributions( probs, rel_ent):
    C = np.zeros(probs.shape)
    for i in range(probs.shape[0]):
        for j in range(probs.shape[1]):
            # the probabilty of a letter times the relative entropy in that sequence column
            C[i][j] = probs[i][j] * rel_ent[j]
    return C
 
 
def get_from_file(fname):
    """ reads data from file and writes values into a nested list"""
    result = []
    with open(fname) as f:
        for line in f:
            innerList = [ch for ch in line.strip()]
            if(innerList):
                result.append(innerList)
    return result


def create_scores( probs , aas):
    """Creates the scores for the plot"""
    scores = []
    for i in range(probs.shape[1]):
        position = []
        for j in range(aas.shape[0]):
            if C[j][i] != 0.0:
                position.append((aas[j],C[j][i]))
        scores.append(position)

    return scores


def sort_scores(scores_to_sort):
    """sorts the scores with regards to the second value in the nested list"""
    sorted_list = []
    for i in scores_to_sort:
        i = sorted(i, key = lambda x: x[1])
        sorted_list.append(i)
    return sorted_list


def show_contributions_table( contributions ):
    """Pretty prints the contributions table with 5 decimals"""
    for i in range(contributions.shape[0]):
        for j in range(contributions.shape[1]):
            print "%.5f" % contributions[i][j],
        print


def most_likely_word(contributions, aas):
    print ''.join(aas[np.argmax(contributions, axis=0)])

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
        font_family = 'DejaVu Sans'

    mpl.rcParams['font.family'] = font_family

    fig, ax = plt.subplots(figsize=(len(all_scores), 5))

    font = FontProperties()
    font.set_size(size)
    font.set_weight('bold')
    
    font.set_family(font_family)

    ax.set_xticks(range(1,len(all_scores)+1))    
    ax.set_yticks(range(0,ymax))
    ax.set_xticklabels(range(1,len(all_scores)+1), rotation=90)
    ax.set_yticklabels(np.arange(0,ymax))
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
                          fontsize=size,
                          color=COLOR_SCHEME[base],
                          ha='center',
                          fontproperties=font,

                         )
            txt.set_path_effects([Scale(0.7, score)])
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
    plt.savefig( filename )

if __name__ == '__main__':


    ## initial setup
    sys.argv += ['demo-input.txt']
    input_file = sys.argv[1]
    plot_file = input_file.split('.')[0] + '.png'
    selo_input = get_from_file(input_file)

    ## convert the input sequences to a numpy array
    SEQs = np.array( selo_input, dtype=str)

    ## converts the array of different letters to a numpy array
    AAs = np.array(different_letters( selo_input ))
    L = len( SEQs[0] ) # the length of each sequence, assumes all the seqs have the same length

    PI = background_freq( AAs, SEQs )
    P = prob_distribution(AAs, SEQs)
    H = entropy(P, PI, SEQs)
    RH = relative_entropy(P, PI, SEQs)
    C = contributions(P, RH)
    #
    print
    print "Contributions:"
    print
    show_contributions_table( C )
    # computing the scores from the contributions and que array of letters
    scores = create_scores(C, AAs)
    sorted_scores = sort_scores(scores)

    print
    print "The most likely word is:"
    print
    most_likely_word(C, AAs)

    print
    print "Storing plot in file: " + plot_file
    print

    ymax = int(ceil(log(max(np.sum(C,axis=0)), 2)))
    draw_logo(sorted_scores, filename=plot_file, ymax=ymax)

    print "BYE"