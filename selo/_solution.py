"""
Used the function from this link as template:

https://gist.github.com/saketkc/45559a011a354a9cca2fbe4e208dde61

"""

#%matplotlib inline
import seaborn
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
from matplotlib import transforms
import matplotlib.patheffects
import numpy as np
import argparse

# Giving the name of the file
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file",
                        required=True,
                        type=str,
                        default=None,
                        help="File containing sequences")
parser.add_argument("--use_back_freq",
                        required=False,
                        action='store_true',
                        default=False,
                        help="If True, uses background frequency")

args = parser.parse_args()


# Opening and reading the file
seqs = []
with open(args.file, 'r') as f:
    for i in f:
        line = i.strip('\n')
        if len(line) != 0:
            seqs.append(line)


# getting all different letters from the sequences
temp = ''
for i in seqs:
    temp += i
temp = ''.join(set(temp))


# Find ing background freq
freq = np.zeros((len(temp)))
for i in range(len(temp)):
    for ii in seqs:
        freq[i] += ii.count(temp[i])
freq = freq / np.sum(freq)
print 'the most probable letter: ', temp[np.argmax(freq)]


# Retreiving color names
colors = []
for name, hex in matplotlib.colors.cnames.iteritems():
    colors.append(name)


# Combining each letter with a color
COLOR_SCHEME = {}
for i in range(len(temp)):
    COLOR_SCHEME[temp[i]] = colors[i]
BASES = list(COLOR_SCHEME.keys())


# making a frequency matrix
m = np.zeros((len(temp), len(seqs[0])))
for i in range(len(seqs[0])):
    for ii in range(len(seqs)):
        m[temp.index(seqs[ii][i])][i] += 1

m1 = m / np.sum(m,axis=0)
m1 = m1.T
m2 = m1 + 0.000000001

print 'the most probable word:'
for i in m1:
    print temp[np.argmax(i)],


# Calculating information content
if args.use_back_freq:
    inf = np.sum(
        m2 * np.log2(m2/np.array(freq)), axis=-1).reshape((-1,1))
else:
    inf = np.log2(len(temp)) + np.sum(m2 * np.log2(m2), axis=-1).reshape((-1,1))
infm = inf * m1



# Making input to the function
inp = []
for i in range(len(seqs[0])):
    inbetween = []
    hk = []
    sorting = []
    sortlet = []
    for ii in range(len(temp)):
        inbetween.append((temp[ii], infm[i][ii]))
    sorting.append([infm[i][ii], temp[ii]])
    sorting.sort()
    for kk in sorting:
        for kkk in inbetween:
            if kk[0] == kkk[1] and kk[1] == kkk[0] and kk != 0:
                hk.append(kkk)
    inp.append(hk)
	






class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def draw_logo(all_scores):
    fig = plt.figure()
    fig.set_size_inches(len(all_scores),5) ###
    ax = fig.add_subplot(111)
    ax.set_xticks(range(len(all_scores)))

    xshift = 0
    trans_offset = transforms.offset_copy(ax.transAxes, 
                                      fig=fig, 
                                      x=0, 
                                      y=0, 
                                      units='points')


    for scores in all_scores:
        yshift = 0

        for base, score in scores:
            txt = ax.text(0, 
                          0, 
                          base, 
                          transform=trans_offset,
                          fontsize=80, 
                          color=COLOR_SCHEME[base],
                          weight='bold',
                          ha='center',
                          family='sans-serif'
                          )
            txt.set_clip_on(False) 
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift = window_ext.height*(score*0.6) ###
            trans_offset = transforms.offset_copy(txt._transform, fig=fig, y=yshift, units='points')
        xshift += window_ext.width*0.74 ###
        trans_offset = transforms.offset_copy(ax.transAxes, fig=fig, x=xshift, units='points')


    ax.set_yticks(range(0,6))###


    seaborn.despine(ax=ax, offset=30, trim=True)
    ax.set_xticklabels(range(1,len(all_scores)+1), rotation=90)
    ax.set_yticklabels(np.arange(0,6,1)) ###
    plt.show()
    #plt.savefig('seqlogo.svg')



draw_logo(inp)


