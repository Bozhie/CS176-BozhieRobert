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

    def getChar(ix, pos):
        """ Returns the character in the 'pos' position of ix suffix. """
        return s[ix + pos]

    def sortSuffixes(i, j, pos):
        """ sorts the suffixes in range i:j in A by the letter in position pos. """

        if (i >= j):
            return

        # tmp_A will be the temporary sorted portion of A for the values from i:j
        tmp_A = np.empty((j - i + 1), dtype=np.int)
        # counter array holding counts for each character in lexicographic order ($ count, A count, ...)
        counter = np.zeros(5, dtype=np.int)

        # C[i] holds the position of the start of a block of each character in A
        C = np.zeros(5, dtype=np.int)
        # tmp_C[i] holds the position of the start of a block of each character in tmp_A
        tmp_C = np.zeros(5, dtype=np.int)


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


    A = np.empty(len(s), dtype=np.int)

    for i in range(len(s)):
        A[i] = i

    sortSuffixes(0, len(s) - 1, 0)

    return list(A)

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
    return list(F)

    pass

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

    pass

def get_occ(L):
    """
    Returns the helper data structure OCC (using the notation from class). OCC should be a dictionary that maps
    string character to a list of integers. If c is a string character and i is an integer, then OCC[c][i] gives
    the number of occurrences of character "c" in the bwt string up to and including index i
    """
    occ = {'$' : [], 'A' : [], 'C' : [], 'G' : [], 'T' : []}
    a,g,t,c,d = 0,0,0,0,0

    for s in L:
        if s == '$':
            d += 1
        elif s == 'A':
            a += 1
        elif s == 'C':
            c += 1
        elif s == 'G':
            g += 1
        else:
            t += 1
        occ['A'].append(a)
        occ['C'].append(c)
        occ['G'].append(g)
        occ['T'].append(t)
        occ['$'].append(d)
    return occ

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

    pos = len(p) - 1
    sp = M[p[pos]]
    if p[pos] == '$':
        ep = M['A'] - 1
    elif p[pos] == 'A':
        ep = M['C'] - 1
    elif p[pos] == 'C':
        ep = M['G'] - 1
    elif p[pos] == 'G':
        ep = M['T'] - 1
    else:
        ep = len(occ['A']) - 1



    while sp <= ep and pos > 0:
        print(sp, ep)
        prev_sp, prev_ep, = sp, ep
        pos -= 1
        c = p[pos]
        sp = M[c] + occ[c][sp - 1]
        ep = M[c] + occ[c][ep] - 1
        if sp > ep:
            return ((prev_sp, prev_ep), len(p) - pos - 1)


    return ((sp, ep), len(p))

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
        
                
        transcriptome = ""
        
        for gene in known_genes:
            for isoform in gene.isoforms:
                iso = ""
                for exon in isoform.exons: 
                    iso = iso + genome_sequence[exon.start:exon.end]
                transcriptome = transcriptome + iso + "$$$$$$"

                
        trans_SA = get_suffix_array(transcriptome)
        
        genome_SA = get_suffix_array(genome_sequence)


        pass


    def simple_align(reads, T):
        #take longest hits, try to extend to hits left and right
        #merge with nearby hits
        #once merged as much as possible, calculate score
        #Take highest scoring match
        A = get_suffix_array(T)
        L = get_bwt(T, A)
        F = get_F(L)
        M = get_M(F)
        occ = get_occ(L)
        alignments = []
        for r in reads:
            rest = r
            potential_aligns = []
            hits = []
            best_so_far = (None, 0)
            maxlen = 0
            while len(rest) > 0:
                result = exact_suffix_matches(rest, M, occ)
                start = result[0][0]
                end = result[0][1]
                if end - start <= 5:
                    for i in range(start, end + 1):
                        hit = (len(rest) - result[1], A[i], result[1])
                        #(position in read, position in genome, length)
                        hits.append(hit)
                rest = rest[:len(rest) - result[1]]

            hits.sort(key=lambda x: x[1])
            i = 0
            while i < len(hits):
                num_extended = 1
                start = hits[i][0]
                genomic_start = hits[i][1]
                mismatches = 0
                while start > 0:
                    start -= 1
                    genomic_start -= 1
                    if r[start] != T[genomic_start]:
                        mismatches += 1
                    if mismatches > 6:
                        break
                if mismatches > 6:
                    i += 1
                    continue
                end = hits[i][0] + hits[i][2] - 1
                curr = i
                next_hit = i + 1
                while next_hit < len(hits) - 1 and hits[next_hit][1] - hits[curr][1] == hits[next_hit][0] - hits[curr][0]:
                    genomic_pos = hits[curr][1] + hits[curr][2]
                    for j in range(1, hits[next_hit][1] - hits[curr][1]):
                        if r[end + j] != T[genomic_pos + j]:
                            mismatches += 1
                    end = hits[next_hit][0] + hits[next_hit][2]
                    curr = next_hit
                    next_hit += 1
                    num_extended += 1
                i += num_extended
                while end < len(r) - 1:
                    end += 1
                    if r[end] != T[genomic_start + end]:
                        mismatches += 1
                if mismatches <= 6:
                    p = ((start, genomic_start, end + 1), mismatches)
                    potential_aligns.append(p)
            best = ()
            min_mismatches = 7
            for x in potential_aligns:
                if x[1] < min_mismatches:
                    best = x
            alignments.append(best)
        return alignments





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
