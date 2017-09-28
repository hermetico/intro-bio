from math import log
import numpy as np
import argparse


def read_sequences_from_file(fname):
    result = []
    with open(fname) as f:
        for line in f:
            innerList =  [ch for ch in line.strip()]
            if(innerList):
                result.append(innerList)
    return result


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
            counter+=1
    return counter


def create_count_matrix(sequences):
    letters = sequences[0].length
    different_letters = different_letters(sequences)
    count_matrix = np.zeros((different_letters, letters))
    for i in range(0, different_letters):
        for j in range(0, letters):
            count_matrix[i][j] = count_apperances(different_letters(j), sequences, i)
    return count_matrix



if __name__ == "__main__":
    # Giving the name of the file
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
                        required=True,
                        type=str,
                        default=None,
                        help="File containing sequences")