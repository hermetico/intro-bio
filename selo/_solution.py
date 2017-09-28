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
        seqs.append(i.strip('\n'))


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
    inf = np.sum(m2 * np.log2(m2/np.array(freq)), axis=-1).reshape((-1,1))
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
    #plt.show()
    plt.savefig('seqlogo.svg')

ALL_SCORES2 = [[('A', 0.01653482213365913),
              ('G', 0.026710097292833978),
              ('C', 0.035613463057111966),
              ('T', 0.057235922770358522)],
             [('C', 0.020055669245080433),
              ('G', 0.023816107228533015),
              ('A', 0.031336983195438178),
              ('T', 0.058913528407423782)],
             [('T', 0.018666958185377256),
              ('G', 0.084001311834197651),
              ('A', 0.093334790926886277),
              ('C', 0.30333807051238043)],
             [('C', 0.0),
              ('G', 0.0),
              ('A', 0.32027512306044359),
              ('T', 0.82203948252180525)],
             [('C', 0.012698627658037786),
              ('A', 0.053334236163758708),
              ('T', 0.096509570201087178),
              ('G', 0.10920819785912497)],
             [('C', 0.0),
              ('G', 0.089472611853783468),
              ('A', 0.1930724782107959),
              ('T', 0.22132698721725386)],
             [('C', 0.020962390607965918),
              ('A', 0.026202988259957396),
              ('G', 0.066380903591892068),
              ('T', 0.07336836712788071)],
             [('G', 0.0),
              ('A', 0.10236420974570831),
              ('C', 0.15354631461856247),
              ('T', 0.29173799777526871)],
             [('G', 0.027681850851852024),
              ('C', 0.089966015268519078),
              ('A', 0.089966015268519078),
              ('T', 0.53287562889815143)],
             [('A', 0.034165612000664765),
              ('C', 0.06833122400132953),
              ('G', 0.072601925501412631),
              ('T', 0.28186629900548432)],
             [('G', 0.0),
              ('A', 0.037325935579058833),
              ('C', 0.23328709736911771),
              ('T', 0.72785574379164719)],
             [('A', 0.017470244196759552),
              ('C', 0.062892879108334396),
              ('G', 0.094339318662501587),
              ('T', 0.19916078384305891)],
             [('G', 0.0),
              ('A', 0.096447131567581681),
              ('C', 0.15844885900388422),
              ('T', 0.48223565783790845)],
             [('G', 0.0),
              ('A', 0.069291952024925829),
              ('C', 0.20787585607477749),
              ('T', 0.46425607856700307)],
             [('G', 0.0),
              ('A', 0.0),
              ('C', 0.21713201856318373),
              ('T', 1.1495224512168551)],
             [('G', 0.0),
              ('A', 0.048934292002649343),
              ('T', 0.27263391258618919),
              ('C', 0.42642740173737281)],
             [('A', 0.0),
              ('G', 0.053607190685875404),
              ('C', 0.2054942309625224),
              ('T', 0.69689347891638032)],
             [('G', 0.0),
              ('A', 0.0),
              ('C', 0.31312908494534769),
              ('T', 0.84220926295645249)],
             [('G', 0.0),
              ('C', 0.068079835765814778),
              ('A', 0.068079835765814778),
              ('T', 1.3207488138568066)],
             [('G', 0.020257705570431345),
              ('A', 0.020257705570431345),
              ('C', 0.048618493369035232),
              ('T', 0.055371061892512348)],
             [('G', 0.0),
              ('A', 0.076286510680262556),
              ('C', 0.20538675952378382),
              ('T', 0.34622339462580698)]]

draw_logo(inp)


