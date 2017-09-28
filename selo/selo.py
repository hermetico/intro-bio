from math import log, ceil

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
        freq[i] = np.sum(seqs == ch)
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


def information_content( probs, aas, rel_ent):
    return log(aas.shape[0]) -rel_ent


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
    """Creates the scores for the plot, and sorts them"""
    scores = []
    for i in range(probs.shape[1]):
        position = []
        for j in range(aas.shape[0]):
            if C[j][i] > 0.0:
                position.append((aas[j],C[j][i]))
        scores.append(sorted(position, key=lambda x:x[1]))

    return scores


def show_table( contributions ):
    """Pretty prints the contributions table with 5 decimals"""
    for i in range(contributions.shape[0]):
        for j in range(contributions.shape[1]):
            print "%.5f" % contributions[i][j],
            #print contributions[i][j],
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
        font_family = 'sans-serif'




    fig, ax = plt.subplots()
    fig.set_size_inches(len(all_scores), 5)  ###
    ax.set_xticks(range(1,len(all_scores)+1))
    ax.set_yticks(range(0, ymax))
    ax.set_xticklabels(range(1,len(all_scores)+1))
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
            yshift = window_ext.height * (score * .6) ###
            trans_offset = transforms.offset_copy(txt._transform, fig=fig,
                                                  y=yshift,
                                                  units='points')
        xshift += window_ext.width * 0.8  ###
        trans_offset = transforms.offset_copy(ax.transAxes, fig=fig,
                                              x=xshift,
                                              units='points')
    plt.savefig( filename )

    #plt.show()


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
    print P
    H = entropy(P, PI, SEQs)
    RH = relative_entropy(P, PI, SEQs)
    C = contributions(P, RH)
    I = information_content(P, AAs, RH)

    # we use the same function, but its not the contributions
    #information_whateva = contributions
    #C = -1 * information_whateva(P, I)
    #
    print
    print "Contributions:"
    print
    show_table( C )
    # computing the scores from the contributions and que array of letters
    scores = create_scores(C, AAs)
    #sorted_scores = sort_scores(scores)

    print
    print "The most likely word is:"
    print
    most_likely_word(C, AAs)

    print
    print "Storing plot in file: " + plot_file
    print

    ymax = int(ceil(max(np.sum(C,axis=0))))
    draw_logo(scores, filename=plot_file, ymax=ymax)
    #draw_logo2(scores)

    print "BYE"