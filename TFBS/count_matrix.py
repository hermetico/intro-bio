from math import log
import numpy as np

def read_sequences_from_file(fname):
    sequences = []
    with open(fname) as f:
        for line in f:
            innerList = line.split()
            if(innerList):
                sequences.append(innerList)
    return sequences

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
    letters = len(sequences[0])
    different_letters1 = different_letters(sequences)
    count_matrix = np.zeros((len(different_letters1), letters))
    for i in range(0, len(different_letters1)):
        for j in range(0, letters):
            count_matrix[i][j] = count_apperances(different_letters(j), sequences, i)
    return count_matrix

readsequences = read_sequences_from_file("input.txt")
result = create_count_matrix(readsequences)
print(result)


             
    
    
    
    
    

