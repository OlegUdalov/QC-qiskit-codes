'''this file contains some usefull functions I made'''
''' OU, 10/2024'''

import sys
import math

#converting int number to a binary represented by a string of ones and zeroes 
def int_2_bin_word(a,n):
    '''converting int number to a binary represented by a string of ones and zeroes 
    converted to int like 5 to a string like '101'. Number of digits in the list is defined by the number n
    Note that the first element of the list represents the highest digit
    
    Parameter:
    __________
    a, int, integer number  to convert
    n, int, number of digits we should have

    Errors:
    __________
    1) if the input number is less than 0, then the -1 is returned and text message is given
    
    Warnings:
    1) if the n is less that number of digit n we set, then the warning is provided.
    
    Returns:
    __________
    string, string of zeroes and ones
    '''
    if a < 0:
        b = -1
        print('int_2_bin_word: Error: number should be positive')
    else:
        b = bin(a)
        s = len(b)
        b = b[slice(2,s)]
        s = s - 2
        if s > n:
            print('int_2_bin_word: Warning: need more digits')
        else:
            if s < n:
                dig = ''
                for i in range(n-s):
                    dig = dig + '0'

                b = dig + b
                    
       
    return b

def plot_counts(counts, bit_size):
    ''' the function plots number of counts for each state after measurements
    The input parameter is the counts object from the result.get_counts(cicr) function. This object is a dictionary with keys being the qubit states and the values corresponding to the key is the number of runs giving such an output. Note that the dictionary is not ordered and the key goes in some random order. 
    The function order the outputs like 0000, 0001, 0010, 0011 etc. The the function just plot the counts as a function of state number.

    Parameter:
    __________
    counts, dictionary with number of counts for each qubits state
    bit_size, int, number of qubits in the circuit

    
    Returns:
    __________
    No returned value, just plot the hystpgram


    '''
    n_hyst = int(math.pow(2,bit_size))
    hyst = []
    for i in range(n_hyst):
        word = int_2_bin_word(i,bit_size)
        if word in counts.keys():
            hyst.append(counts[word])
        else:
            hyst.append(0)
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    xpoints = range(n_hyst)
    
    plt.plot(xpoints, hyst)
    plt.show()

def plot_counts_2(counts_1, counts_2, bit_size):
    ''' the function plots number of counts for each state after measurements
    IN contrast to previous function it plots 2 graph on the same plot for comparison
    The input parameter is the counts object from the result.get_counts(cicr) function. This object is a dictionary with keys being the qubit states and the values corresponding to the key is the number of runs giving such an output. Note that the dictionary is not ordered and the key goes in some random order. 
    The function order the outputs like 0000, 0001, 0010, 0011 etc. The the function just plot the counts as a function of state number.

    Parameter:
    __________
    counts, dictionary with number of counts for each qubits state
    bit_size, int, number of qubits in the circuit

    
    Returns:
    __________
    No returned value, just plot the hystpgram


    '''
    n_hyst = int(math.pow(2,bit_size))
    hyst_1 = []
    for i in range(n_hyst):
        word = int_2_bin_word(i,bit_size)
        if word in counts_1.keys():
            hyst_1.append(counts_1[word])
        else:
            hyst_1.append(0)

    hyst_2 = []
    for i in range(n_hyst):
        word = int_2_bin_word(i,bit_size)
        if word in counts_2.keys():
            hyst_2.append(counts_2[word])
        else:
            hyst_2.append(0)
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    xpoints = range(n_hyst)
    
    plt.plot(xpoints, hyst_1, label = ['1'])
    plt.plot(xpoints, hyst_2, label = ['2'])
    plt.legend()
    plt.show()
        
        