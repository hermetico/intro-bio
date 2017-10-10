import numpy as np
import argparse
import seq_logo_plot

def read_sequences_from_file( fname ):
    """Reads sequences from fasta file"""
    sequences = []
    with open(fname) as f:
        for line in f:
            # skips lines with sequence name
            if line[0] == '>':
                continue

            sequence = list(line.rstrip().lower())
            # skips empty lines
            if(sequence):
                sequences.append(sequence)
    return np.array(sequences)





class EMA(object):

    def __init__(self, setup):
        self.setup = setup

        self.SEQs = setup['seqs']
        self.len_motifs = setup['lengths']
        self.pseudo_counts = setup.get('pseudo_counts', None)
        self.threshold = setup.get('threshold', 0.01)
        self.AAs = np.array(['a', 'c', 'g', 't'])
        self.P0 = self.__compute_background_freq__(self.SEQs, self.AAs)
        #dict with the index associated to each aas
        self.aas_index = { a:i for i, a in enumerate(self.AAs)}


    def __compute_background_freq__(self,  seqs, aas ):
        """Computes the background frequency of the given aas in a list of sequences
        returns a 2D array
        """
        totals = np.zeros(aas.shape[0])
        for i, ch in enumerate(aas):
            totals[i] = (seqs == ch).sum()
        return totals[:, np.newaxis] / totals.sum()


    def __add_initial_positions__(self, z, guesses=None):
        if guesses is None:
            # pics random integers from 0 to len((sequence - len_motif) + 1)
            # as many times as sequences
            #guesses = np.random.choice(z.shape[1], z.shape[0])
            # try with the last one
            guesses = [ z.shape[1] -1 for x in range(z.shape[0])]
            # try with the last first one
            #guesses = [ 0 for x in range(z.shape[0])]

        for i, position in enumerate(guesses):
                z[i][position] = 1
        return z


    def __pwm_count__(self, sequences, aas):
        """Computes the pwm count matrix"""
        count_matrix = np.zeros((aas.shape[0], sequences.shape[1]))
        for i, ch in enumerate(aas):
            count_matrix[i, :] = (sequences == ch).sum(axis=0)

        return count_matrix


    def __pwm_freq__(self, sequences, aas):
        """Computes the pwm freq with or without pseudocounts"""
        counts = self.__pwm_count__(sequences, aas)
        if self.pseudo_counts is not None:
            # adding the pseudocounts to the background frequencies
            back_freq_column = self.P0 * self.pseudo_counts
            # increases each cell of freqs by the value corresponding to its back_freq_column row
            counts +=  back_freq_column

        # returns everything normalized
        return counts / counts.sum(axis=0)


    def __extract_motifs__(self, seqs, z, uniques=True):
        """Extracts unique motifs, using the probability for a certain first position"""
        motifs = []
        if uniques:
            starting = enumerate(z.argmax(axis=1)) # the maximum sarting
        else:
            starting = np.argwhere(z > 0)

        # creates the candidate motif
        for i, start in starting:
            motifs.append(seqs[i, start:start + self.len_motifs])

        return np.array(motifs)


    def __compute_tht__(self, seqs, z, aas):
        """
        """
        starting = np.argwhere(z > 0)  # all above 0

        tht = np.zeros((aas.shape[0], self.len_motifs))
        for i, start in starting:
            motif = seqs[i, start:start + self.len_motifs]
            prob = z[i][start]
            for j, letter in enumerate(motif):
                ## adds to the correct row (for a particular letter)
                ## the value that motif starting on that z position
                tht[ self.aas_index[letter]][j] += prob
        return tht / tht.sum(axis=0)


    def __meme_extension23__(self, seqs):
        """Looking for a good starting point (slide 54)"""

        ## assume all the possible starting points
        z = np.ones((seqs.shape[0],
                     (seqs.shape[1] - self.len_motifs) + 1))

        ## get all the motifs
        motifs = self.__extract_motifs__(self.SEQs, z, uniques=False)
        th = self.__pwm_freq__(motifs, self.AAs)

        ## run one EM step
        z, _ = self.__iterate__(self.SEQs, z, th)

        ## chose the model with highest likelihood
        # new z
        nz = np.zeros(z.shape)

        return self.__add_initial_positions__(nz, z.argmax(axis=1))


    def compute_ema(self, setup=None):
        """Computes the expectation maximization algorithm"""
        if setup is not None:
            self.__init__(setup)

        self.Z = np.zeros((self.SEQs.shape[0],
           (self.SEQs.shape[1] - self.len_motifs) + 1))

        # initial guess
        self.Z = self.__add_initial_positions__(self.Z)
        #self.P0 = self.__compute_background_freq__(self.SEQs[:,: self.Z.shape[1] - 1],self.AAs)
        #self.Z = self.__meme_extension23__(self.SEQs)
        # first motif representation
        motifs = self.__extract_motifs__(self.SEQs, self.Z)
        self.Th = self.__pwm_freq__(motifs, self.AAs)


        diff = np.inf
        iters = 0

        while diff > 0:#self.threshold: ## threshold is important!!
            # compute new Z
            Zt, Tht = self.__iterate__(self.SEQs, self.Z, self.Th)

            # compute difference
            diff = np.sum(np.absolute(self.Th - Tht))
            #diff2 = np.sum(np.absolute(self.Z - Zt))

            # update Z-1 and Th-1
            self.Z = Zt
            self.Th = Tht

            iters += 1
            #print "diff ", diff
            #print "diff2 ", diff2
        # updates the list of motifs accordingly
        self.motifs = self.__extract_motifs__(self.SEQs, self.Z)
        print iters


    def most_likely_sequence(self):
        """Returns the most likely sequence"""
        return ''.join(self.AAs[self.Th.argmax( axis=0)])

    def __iterate__(self, seqs, z, th):

        # Z on iteration t
        zt = np.zeros(z.shape)

        # E-STEP
        # each sequence
        for i in range(z.shape[0]): # sequence
            for j in range(z.shape[1]): # position

                p = 1L
                ## before background
                # from 0 to j
                for k in range(j):
                    Cik = seqs[i][k]
                    # product with the background of character Cik
                    p *= self.P0[self.aas_index[Cik]]

                ## foreground
                # from the beginning of the motifs to its end
                for k in range(j, j + self.len_motifs):
                    Cik = seqs[i][k]
                    # product with the foreground of character Cik
                    p *= th[self.aas_index[Cik]][k - j]

                ## after background
                # from end of motifs to end of sequence
                for k in range(j + self.len_motifs, seqs.shape[1]):
                    Cik = seqs[i][k]
                    # product with the background of character Cik
                    p *= self.P0[self.aas_index[Cik]]

                # update Z of
                zt[i][j] = p

        # normalize all Z row by row
        zt = zt / zt.sum(axis=1)[:, np.newaxis]

        ### The bug is in here
        # M-STEP update the motif model
        ##motifs = self.__extract_motifs__(seqs, zt)
        # theta on iteration t
        ###tht = self.__pwm_freq__(motifs, self.AAs)


        ## M-STEP
        # updating th accordingly
        tht = self.__compute_tht__(seqs, zt, self.AAs)

        # returns both zt and tht
        return zt, tht




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
                        required=True,
                        type=str,
                        default=None,
                        help="Fasta file containing sequences")

    parser.add_argument("-l", "--lengths",
                        required=False,
                        type=int,
                        default=19,
                        help="Length of most likely binding sequences base-pairs")

    parser.add_argument("-th", "--threshold",
                        required=False,
                        type=int,
                        default=None,
                        help="Threshold to end iterations")

    parser.add_argument("-o", "--output",
                        required=False,
                        type=str,
                        default='seqlogo.png',
                        help="Output sequence logo filename")

    args = parser.parse_args()

    SEQS = read_sequences_from_file(args.file)

    # setup
    setup = {}
    setup['lengths'] = args.lengths
    setup['seqs'] = SEQS

    if args.threshold is not None:
        setup['threshold'] = args.threshold


    ploter = seq_logo_plot.Logo_Ploter()
    ema = EMA(setup)
    ema.compute_ema()

    # motifs
    print "Motifs:"
    for m in ema.motifs:
        print ''.join(m)

    print
    print "Consensus sequence:"
    print ema.most_likely_sequence()

    pwm = ema.Th
    ## uses relative entropy
    ploter.compute_logo(pwm, ema.motifs, ema.P0, ema.AAs, method=ploter.relative_entropy, filename=args.output)
    ## uses normal entropy
    #ploter.compute_logo(pwm, ema.motifs, ema.P0, ema.AAs, filename=args.output)

