import numpy as np
import argparse
import seq_logo_plot
import math
FREQs = {'a': .2, 't': .2, 'c': .3, 'g': .3}


def read_sequences_from_file(fname):
    sequences = []
    with open(fname) as f:
        for line in f:
            innerList = [ ch for ch in line.rstrip()]
            if(innerList):
                sequences.append(innerList)
    return np.array(sequences)


def create_pwm_count_matrix(sequences, aas):
    """Computes the pwm count matrix"""
    count_matrix = np.zeros((aas.shape[0], sequences.shape[1]))
    for i, ch in enumerate(aas):
        count_matrix[i,:] = (sequences == ch).sum(axis=0)

    return count_matrix


def create_pwm_freq_matrix( pwm_count, pseudo_counts=None, back_freq=None):
    """Computes the pwm freq with or without pseudocounts"""
    if pseudo_counts is not None and back_freq is not None:
        # adding the pseudocounts to the background frequencies
        back_freq_column = back_freq * pseudo_counts
        # increases each cell of freqs by the value corresponding to its back_freq_column row
        freqs = pwm_count + back_freq_column
        # returns everything normalized
        return freqs / freqs.sum(axis=0)

    # returns everything normalized
    return pwm_count / pwm_count.sum(axis=0)



def create_relative_entropy(probs, back_probs):
    """Computes the relative entropy"""
    if np.any(probs == 0):
        RH = np.zeros(probs.shape[1])
        for j in range(probs.shape[1]):
            accum = 0.0
            for i in range(back_probs.shape[0]):
                if probs[i][j] != 0.0 and back_probs[i] != 0.0:  # avoiding zeros
                    # the sum in a sequence column of the probability of each letter times the logarithm in base 2
                    # of the probability of a letter divided by its background freq
                    accum += probs[i][j] * np.log2(probs[i][j] / back_probs[i])
            RH[j] = accum
        return RH

    return np.sum( probs *  np.log2( probs / back_probs),
        axis=0).ravel()



def print_matrix( m, precision=2):
    """Prints the matrix with the selected precision"""
    for i in range(m.shape[0]):
        print ", ".join(["%.*f" %(precision, num) for num in m[i]])
    print("")



if __name__ == "__main__":

    ploter = seq_logo_plot.Logo_Ploter()

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
                        required=True,
                        type=str,
                        default=None,
                        help="File containing sequences")

    args = parser.parse_args()

    sequences = read_sequences_from_file(args.file)

    AAs = np.unique(sequences)
    PI = np.array([FREQs[a] for a in AAs]) # one array
    PI = PI[:, np.newaxis] #2D array with one column
    PWM_count = create_pwm_count_matrix(sequences, AAs)

    print("PI:")
    for i in range(PI.shape[0]):
        print("%s: %f"% (AAs[i], PI[i][0]))

    print("\nPWM count")
    print(PWM_count)
    print("\nPWM freq")


    for pseudocount in [0, 1, 3, 5]:
        print("With pseudocounts %i" %(pseudocount))
        if pseudocount == 0:
            pwm = create_pwm_freq_matrix(PWM_count)
        else:
            pwm = create_pwm_freq_matrix(PWM_count, pseudo_counts=pseudocount, back_freq=PI)

        assert(np.sum(pwm.sum(axis=0)) == pwm.shape[1]) # checks everything is normalized
        print_matrix(pwm)
        """
        RH = create_relative_entropy(pwm, PI)
        C = pwm * RH  # contributions

        scores = ploter.create_scores(C, AAs)

        ploter.draw_logo(C, AAs, filename="seq-logo%i.png" %(pseudocount))
        """


    probs = create_pwm_freq_matrix(PWM_count)
    RH = create_relative_entropy(probs, PI)
    C = probs * RH  # contributions

    scores = ploter.create_scores(C, AAs)

    ploter.draw_logo(C, AAs, filename="seq-logo.png")





