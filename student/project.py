""" 
    RNA Alignment Assignment
    
    Implement each of the functions below using the algorithms covered in class.
    You can construct additional functions and data structures but you should not
    change the functions' APIs.

    You will be graded on the helper function implementations as well as the RNA alignment, although
    you do not have to use your helper function.
    
    *** Make sure to comment out any print statement so as not to interfere with the grading script
"""

import sys # DO NOT EDIT THIS
from shared import *
import numpy

ALPHABET = [TERMINATOR] + BASES

def get_suffix_array(s):
    """
    Naive implementation of suffix array generation (0-indexed). You do not have to implement the
    KS Algorithm. Make this code fast enough so you have enough time in Aligner.__init__ (see bottom).

    Input:
        s: a string of the alphabet ['A', 'C', 'G', 'T'] already terminated by a unique delimiter '$'

    Output: list of indices representing the suffix array

    >>> get_suffix_array('GATAGACA$')
    [8, 7, 5, 3, 1, 6, 4, 0, 2]
    """

    def getChar(ix, pos):
        """ Returns the character in the 'pos' position of ix suffix. """
        return s[ix + pos]

    def sortSuffixes(i, j, pos):
        """ sorts the suffixes in range i:j in A by the letter in position pos. """

        if (i >= j):
            return

        # tmp_A will be the temporary sorted portion of A for the values from i:j
        tmp_A = numpy.empty((j - i + 1), dtype=numpy.int)
        # counter array holding counts for each character in lexicographic order ($ count, A count, ...)
        counter = numpy.zeros(5, dtype=numpy.int)

        # C[i] holds the position of the start of a block of each character in A
        C = numpy.zeros(5, dtype=numpy.int)
        # tmp_C[i] holds the position of the start of a block of each character in tmp_A
        tmp_C = numpy.zeros(5, dtype=numpy.int)


        for k in range(i, j + 1):
            if (getChar(A[k], pos) == "$"):
                counter[0] += 1
            elif (getChar(A[k], pos) == 'A'):
                counter[1] += 1
            elif (getChar(A[k], pos) == 'C'):
                counter[2] += 1
            elif (getChar(A[k], pos) == 'G'):
                counter[3] += 1
            elif (getChar(A[k], pos) == 'T'):
                counter[4] += 1


        # Fill out C and tmp_C using the counts
        for k in range(0,len(counter)):
            if k == 0:
                C[k] = i
                tmp_C[k] = 0
            else:
                C[k] = C[k-1] + counter[k-1]
                tmp_C[k] = tmp_C[k-1] + counter[k-1]


        for k in range(i, j+1):
            if (getChar(A[k], pos) == '$'):
                tmp_A[tmp_C[0]] = A[k]
                tmp_C[0] += 1
            elif (getChar(A[k], pos) == 'A'):
                tmp_A[tmp_C[1]] = A[k]
                tmp_C[1] += 1
            elif (getChar(A[k], pos) == 'C'):
                tmp_A[tmp_C[2]] = A[k]
                tmp_C[2] += 1
            elif (getChar(A[k], pos) == 'G'):
                tmp_A[tmp_C[3]] = A[k]
                tmp_C[3] += 1
            elif (getChar(A[k], pos) == 'T'):
                tmp_A[tmp_C[4]] = A[k]
                tmp_C[4] += 1

        A[i:j + 1] = tmp_A[0:len(tmp_A) + 1]

        new_pos = pos + 1
        for k in range(0, len(C)):
            if (counter[k] > 1):
                if k == len(C)-1:
                    sortSuffixes(C[k], j, new_pos)
                else:
                    sortSuffixes(C[k], C[k+1] - 1, new_pos)


    A = numpy.empty(len(s), dtype=numpy.int)

    for i in range(len(s)):
        A[i] = i

    sortSuffixes(0, len(s) - 1, 0)

    return A
    
    pass

def get_bwt(s, sa):
    """
    Input:
        s: a string terminated by a unique delimiter '$'
        sa: the suffix array of s

    Output:
        L: BWT of s as a string
    """
    pass

def get_F(L):
    """
    Input: L = get_bwt(s)

    Output: F, first column in Pi_sorted
    """
    pass

def get_M(F):
    """
    Returns the helper data structure M (using the notation from class). M is a dictionary that maps character
    strings to start indices. i.e. M[c] is the first occurrence of "c" in F.

    If a character "c" does not exist in F, you may set M[c] = -1
    """
    pass

def get_occ(L):
    """
    Returns the helper data structure OCC (using the notation from class). OCC should be a dictionary that maps 
    string character to a list of integers. If c is a string character and i is an integer, then OCC[c][i] gives
    the number of occurrences of character "c" in the bwt string up to and including index i
    """
    pass

def exact_suffix_matches(p, M, occ):
    """
    Find the positions within the suffix array sa of the longest possible suffix of p 
    that is a substring of s (the original string).
    
    Note that such positions must be consecutive, so we want the range of positions.

    Input:
        p: the pattern string
        M, occ: buckets and repeats information used by sp, ep

    Output: a tuple (range, length)
        range: a tuple (start inclusive, end exclusive) of the indices in sa that contains
            the longest suffix of p as a prefix. range=None if no indices matches any suffix of p
        length: length of the longest suffix of p found in s. length=0 if no indices matches any suffix of p

        An example return value would be ((2, 5), 7). This means that p[len(p) - 7 : len(p)] is
        found in s and matches positions 2, 3, and 4 in the suffix array.

    >>> s = 'ACGT' * 10 + '$'
    >>> sa = get_suffix_array(s)
    >>> sa
    [40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0, 37, 33, 29, 25, 21, 17, 13, 9, 5, 1, 38, 34, 30, 26, 22, 18, 14, 10, 6, 2, 39, 35, 31, 27, 23, 19, 15, 11, 7, 3]
    >>> L = get_bwt(s, sa)
    >>> L
    'TTTTTTTTTT$AAAAAAAAAACCCCCCCCCCGGGGGGGGGG'
    >>> F = get_F(L)
    >>> F
    '$AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT'
    >>> M = get_M(F)
    >>> sorted(M.items())
    [('$', 0), ('A', 1), ('C', 11), ('G', 21), ('T', 31)]
    >>> occ = get_occ(L)
    >>> type(occ) == dict, type(occ['$']) == list, type(occ['$'][0]) == int
    (True, True, True)
    >>> occ['$']
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    >>> exact_suffix_matches('ACTGA', M, occ)
    ((1, 11), 1)
    >>> exact_suffix_matches('$', M, occ)
    ((0, 1), 1)
    >>> exact_suffix_matches('AA', M, occ)
    ((1, 11), 1)
    """
    pass

MIN_INTRON_SIZE = 20
MAX_INTRON_SIZE = 10000

class Aligner:
    def __init__(self, genome_sequence, known_genes):
        """
        Initializes the aligner. Do all time intensive set up here. i.e. build suffix array.

        genome_sequence: a string (NOT TERMINATED BY '$') representing the bases of the of the genome
        known_genes: a python set of Gene objects (see shared.py) that represent known genes. You can get the isoforms 
                     and exons from a Gene object

        Time limit: 500 seconds maximum on the provided data. Note that our server is probably faster than your machine, 
                    so don't stress if you are close. Server is 1.25 times faster than the i7 CPU on my computer

        """

        pass

    def align(self, read_sequence):
        """
        Returns an alignment to the genome sequence. An alignment is a list of pieces. 
        Each piece consists of a start index in the read, a start index in the genome, and a length 
        indicating how many bases are aligned in this piece. Note that mismatches are count as "aligned".

        Note that <read_start_2> >= <read_start_1> + <length_1>. If your algorithm produces an alignment that 
        violates this, we will remove pieces from your alignment arbitrarily until consecutive pieces 
        satisfy <read_start_2> >= <read_start_1> + <length_1>

        Return value must be in the form (also see the project pdf):
        [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        If no good matches are found: return the best match you can find or return []

        Time limit: 0.5 seconds per read on average on the provided data.
        """
        pass
