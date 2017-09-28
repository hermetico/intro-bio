from math import log
import numpy as np
import argparse



FREQs = {'a': .2, 't': 0.2, 'c': .3, 'g': .3}

def read_sequences_from_file(fname):
    sequences = []
    with open(fname) as f:
        for line in f:
            innerList = [ ch for ch in line.rstrip()]
            if(innerList):
                sequences.append(innerList)
    return np.array(sequences)


def different_letters( lines ):
    seen = []
    for line in lines:
        for ch in line:
            if ch not in seen:
                seen.append(ch)
    return seen


def count_apperances(letter, sequences, position):
    counter = 0
    for sequence in sequences():
        if(sequence[position] == letter):
            counter += 1
    return counter


def create_pwm_count_matrix(sequences, aas):
    count_matrix = np.zeros((aas.shape[0], sequences.shape[1]))
    for i, ch in enumerate(aas):
        count_matrix[i,:] = (sequences == ch).sum(axis=0)

    return np.matrix(count_matrix)

def create_pwm_freq_matrix( pwm_count, pseudo_count=None, back_freq=None):
    if pseudo_count is not None and back_freq is not None:
        # creates a 2D array with one column that contains the data and multiplies its rows by pseudo_count value
        back_freq_column = back_freq * pseudo_count
        # increases each cell of freqs by the value corresponding to its back_freq_column row
        freqs = pwm_count + back_freq_column
        return freqs / freqs.sum(axis=0)

    return pwm_count / pwm_count.sum(axis=0)

def print_matrix( m, precision=2):
    for i in range(m.shape[0]):
        print ", ".join(["%.*f" %(precision, num) for num in m[i].getA1()])
    print("")



if __name__ == "__main__":
    # Giving the name of the file
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
                        required=True,
                        type=str,
                        default=None,
                        help="File containing sequences")

    sequences = read_sequences_from_file("input.txt")

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
    print("With pseudocounts 0")
    pwm = create_pwm_freq_matrix(PWM_count)
    assert(np.sum(pwm.sum(axis=0)) == pwm.shape[1]) # checks everything is normalized
    print_matrix(pwm)
    print("With pseudocounts 1")
    pwm = create_pwm_freq_matrix(PWM_count, pseudo_count=1, back_freq=PI)
    assert (np.sum(pwm.sum(axis=0)) == pwm.shape[1])
    print_matrix(pwm)
    print("With pseudocounts 2")
    pwm = create_pwm_freq_matrix(PWM_count, pseudo_count=2, back_freq=PI)
    assert (np.sum(pwm.sum(axis=0)) == pwm.shape[1])
    print_matrix(pwm)

    print("With pseudocounts 4")
    pwm = create_pwm_freq_matrix(PWM_count, pseudo_count=4, back_freq=PI)
    assert (np.sum(pwm.sum(axis=0)) == pwm.shape[1])
    print_matrix(pwm)

