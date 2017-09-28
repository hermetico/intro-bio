from math import log
import numpy as np

def read_sequences_from_file:
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

