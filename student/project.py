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
import numpy as np

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
    pass

def get_bwt(s, sa):
    
    L = ""
    for i in range(len(sa)):
        L += s[sa[i] - 1 % len(s)]
    return L

def get_F(L):
    
    a_count, c_count, g_count, t_count, dollar = 0,0,0,0,0
    for s in L:
        if s == 'A':
            a_count += 1
        elif s == 'C':
            c_count += 1
        elif s == 'G':
            g_count += 1
        elif s == 'T':
            t_count += 1
        else:
            dollar += 1
    F = np.empty(len(L), str)
    i = 0
    for j in range(dollar):
        F[i] = '$'
        i += 1
    for j in range(a_count):
        F[i] = 'A'
        i += 1
    for j in range(c_count):
        F[i] = 'C'
        i += 1
    for j in range(g_count):
        F[i] = 'G'
        i += 1
    for j in range(t_count):
        F[i] = 'T'
        i += 1
    return F
        
def get_M(F):
    
    def bin_search(c):
        lo = 0
        hi = len(F)
        while lo < hi - 1:
            mid = (lo + hi)//2 
            if (F[mid] < c):
                lo = mid
            else:
                hi = mid
        return lo + 1
    return {'$' : 0, 'A' : 1, 'C' : bin_search('C'), 'G' : bin_search('G'), 'T': bin_search('T')}


def get_occ(L):
    """
    Returns the helper data structure OCC (using the notation from class). OCC should be a dictionary that maps 
    string character to a list of integers. If c is a string character and i is an integer, then OCC[c][i] gives
    the number of occurrences of character "c" in the bwt string up to and including index i
    """
    occ = {'$' : [], 'A' : [], 'C' : [], 'G' : [], 'T' : []}
    a,g,t,c,d = 0,0,0,0,0
    for c in L:
        if c == '$':
            d = 1
        elif c == 'A':
            d = 1
        elif c == 'C':
            d = 1
        elif c == 'G':
            d = 1
        else c == 'T':
            d = 1
        occ['A'].append(occ[len(occ['A']) - 1] + a)
        occ['C'].append(occ[len(occ['C']) - 1] + c)
        occ['G'].append(occ[len(occ['G']) - 1] + g)
        occ['T'].append(occ[len(occ['t']) - 1] + t)
        occ['$'].append(occ[len(occ['$']) - 1] + d)
    return occ

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
    
    pos = len(p) - 1
    sp = M[p[pos]]
    if p[pos] == '$':
        ep = M['A'] - 1
    elif sp == 'A':
        ep = M['C'] - 1
    elif sp == 'C':
        ep = M['G'] - 1
    elif sp == 'G':
        ep = M['T'] - 1
    else 
        ep = len(occ['A']) - 1
        
    while sp <= ep and pos > 0:
        pos -= 1
        c = p[pos]
        sp = M[c] + occ[c][sp - 1]
        ep = M[c] + occ[c][ep] - 1
    
    return ((sp, ep), len(p) - pos))
        

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
