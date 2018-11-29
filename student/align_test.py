
# coding: utf-8

# In[62]:


import numpy as np
from math import ceil, log2
from project import *
from shared import *


# In[80]:


gene1 = "GCTCTCTGCTCCTCCTGTTCGACAGTCAGCCGCATCTTCTTTTGCGTCGCCAGGTGAAGACGGGCGGAGAGAAACCCGGGAGGCTAGGGACGGCCTGAAGGCGGCAGGGGCGGGCGCAGGCCGGATGTGTTCGCGCCGCTGCGGGGTGGGCCCGGGCGGCCTCCGCATTGCAGGGGCGGGCGGAGGACGTGATGCGGCGCGGGCTGGGCATGGAGGCCTGGTGGGGGAGGGGAGGGGAGGCGTGTGTGTCGGCCGGGGCCACTAGGCGCTCACTGTTCTCTCCCTCCGCGCAGCCGAGCCACATCGCTCAGACACCATGGGGAAGGTGAAGGTCGGAGTCAACGGGTGAGTTCGCGGGTGGCTGGGGGGCCCTGGGCTGCGACCGCCCCCGAACCGCGTCTACGAGCCTTGCGGGCTCCGGGTCTTTGCAGTCGTATGGGGGCAGGGTAGCTGTTCCCCGCAAGGAGAGCTCAAGGTCAGCGCTCGGACCTGGCGGAGCCCCGCACCCAGGCTGTGGCGCCCTGTGCAGCTCCGCCCTTGCGGCGCCATCTGCCCGGAGCCTCCTTCCCCTAGTCCCCAGAAACAGGAGGTCCCTACTCCCGCCCGAGATCCCGACCCGGACCCCTAGGTGGGGGACGCTTTCTTTCCTTTCGCGCTCTGCGGGGTCACGTGTCGCAGAGGAGCCCCTCCCCCACGGCCTCCGGCACCGCAGGCCCCGGGATGCTAGTGCGCAGCGGGTGCATCCCTGTCCGGATGCTGCGCCTGCGGTAGAGCGGCCGCCATGTTGCAACCGGGAAGGAAATGAATGGGCAGCCGTTAGGAAAGCCTGCCGGTGACTAACCCTGCGCTCCTGCCTCGATGGGTGGAGTCGCGTGTGGCGGGGAAGTCAGGTGGAGCGAGGCTAGCTGGCCCGATTTCTCCTCCGGGTGATGCTTTTCCTAGATTATTCTCTGGTAAATCAAAGAAGTGGGTTTATGGAGGTCCTCTTGTGTCCCCTCCCCGCAGAGGTGTGGTGGCTGTGGCATGGTGCCAAGCCGGGAGAAGCTGAGTCATGGGTAGTTGGAAAAGGACATTTCCACCGCAAAATGGCCCCTCTGGTGGTGGCCCCTTCCTGCAGCGCCGGCTCACCTCACGGCCCCGCCCTTCCCCTGCCAGCCTAGCGTTGACCCGACCCCAAAGGCCAGGCTGTAAATGTCACCGGGAGGATTGGGTGTCTGGGCGCCTCGGGGAACCTGCCCTTCTCCCCATTCCGTCTTCCGGAAACCAGATCTCCCACCGCACCCTGGTCTGAGGTTAAATATAGCTGCTGACCTTTCTGTAGCTGGGGGCCTGGGCTGGGGCTCTCTCCCATCCCTTCTCCCCACACACATGCACTTACCTGTGCTCCCACTCCTGATTTCTGGAAAAGAGCTAGGAAGGACAGGCAACTTGGCAAATCAAAGCCCTGGGACTAGGGGGTTAAAATACAGCTTCCCCTCTTCCCACCCGCCCCAGTCTCTGTCCCTTTTGTAGGAGGGACTTAGAGAAGGGGTGGGCTTGCCCTGTCCAGTTAATTTCTGACCTTTACTCCTGCCCTTTGAGTTTGATGATGCTGAGTGTACAAGCGTTTTCTCCCTAAAGGGTGCAGCTGAGCTAGGCAGCAGCAAGCATTCCTGGGGTGGCATAGTGGGGTGGTGAATACCATGTACAAAGCTTGTGCCCAGACTGTGGGTGGCAGTGCCCCACATGGCCGCTTCTCCTGGAAGGGCTTCGTATGACTGGGGGTGTTGGGCAGCCCTGGAGCCTTCAGTTGCAGCCATGCCTTAAGCCAGGCCAGCCTGGCAGGGAAGCTCAAGGGAGATAAAATTCAACCTCTTGGGCCCTCCTGGGGGTAAGGAGATGCTGCATTCGCCCTCTTAATGGGGAGGTGGCCTAGGGCTGCTCACATATTCTGGAGGAGCCTCCCCTCCTCATGCCTTCTTGCCTCTTGTCTCTTAGATTTGGTCGTATTGGGCGCCTGGTCACCAGGGCTGCTTTTAACTCTGGTAAAGTGGATATTGTTGCCATCAATGACCCCTTCATTGACCTCAACTACATGGTGAGTGCTACATGGTGAGCCCCAAAGCTGGTGTGGGAGGAGCCACCTGGCTGATGGGCAGCCCCTTCATACCCTCACGTATTCCCCCAGGTTTACATGTTCCAATATGATTCCACCCATGGCAAATTCCATGGCACCGTCAAGGCTGAGAACGGGAAGCTTGTCATCAATGGAAATCCCATCACCATCTTCCAGGAGTGAGTGGAAGACAGAATGGAAGAAATGTGCTTTGGGGAGGCAACTAGGATGGTGTGGCTCCCTTGGGTATATGGTAACCTTGTGTCCCTCAATATGGTCCTGTCCCCATCTCCCCCCCACCCCCATAGGCGAGATCCCTCCAAAATCAAGTGGGGCGATGCTGGCGCTGAGTACGTCGTGGAGTCCACTGGCGTCTTCACCACCATGGAGAAGGCTGGGGTGAGTGCAGGAGGGCCCGCGGGAGGGGAAGCTGACTCAGCCCTGCAAAGGCAGGACCCGGGTTCATAACTGTCTGCTTCTCTGCTGTAGGCTCATTTGCAGGGGGGAGCCAAAAGGGTCATCATCTCTGCCCCCTCTGCTGATGCCCCCATGTTCGTCATGGGTGTGAACCATGAGAAGTATGACAACAGCCTCAAGATCATCAGGTGAGGAAGGCAGGGCCCGTGGAGAAGCGGCCAGCCTGGCACCCTATGGACACGCTCCCCTGACTTGCGCCCCGCTCCCTCTTTCTTTGCAGCAATGCCTCCTGCACCACCAACTGCTTAGCACCCCTGGCCAAGGTCATCCATGACAACTTTGGTATCGTGGAAGGACTCATGGTATGAGAGCTGGGGAATGGGACTGAGGCTCCCACCTTTCTCATCCAAGACTGGCTCCTCCCTGCCGGGGCTGCGTGCAACCCTGGGGTTGGGGGTTCTGGGGACTGGCTTTCCCATAATTTCCTTTCAAGGTGGGGAGGGAGGTAGAGGGGTGATGTGGGGAGTACGCTGCAGGGCCTCACTCCTTTTGCAGACCACAGTCCATGCCATCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGCGGGGCTCTCCAGAACATCATCCCTGCCTCTACTGGCGCTGCCAAGGCTGTGGGCAAGGTCATCCCTGAGCTGAACGGGAAGCTCACTGGCATGGCCTTCCGTGTCCCCACTGCCAACGTGTCAGTGGTGGACCTGACCTGCCGTCTAGAAAAACCTGCCAAATATGATGACATCAAGAAGGTGGTGAAGCAGGCGTCGGAGGGCCCCCTCAAGGGCATCCTGGGCTACACTGAGCACCAGGTGGTCTCCTCTGACTTCAACAGCGACACCCACTCCTCCACCTTTGACGCTGGGGCTGGCATTGCCCTCAACGACCACTTTGTCAAGCTCATTTCCTGGTATGTGGCTGGGGCCAGAGACTGGCTCTTAAAAAGTGCAGGGTCTGGCGCCCTCTGGTGGCTGGCTCAGAAAAAGGGCCCTGACAACTCTTTTCATCTTCTAGGTATGACAACGAATTTGGCTACAGCAACAGGGTGGTGGACCTCATGGCCCACATGGCCTCCAAGGAGTAAGACCCCTGGACCACCAGCCCCAGCAAGAGCACAAGAGGAAGAGAGAGACCCTCACTGCTGGGGAGTCCCTGCCACACTCAGTCCCCCACCACACTGAATCTCCCCTCCTCACAGTTGCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGCCGCACCTTGTCATGTACCATCAATAAAGTACCCTGTGCTCAACCAGTTA$"
#gene = "ATCGCGCG$"
#gene = "ACAAGTACGTGGTGCTGGGCTATAAAAACAAATGGAACTGCGCTGCAGTTGCTTTAGTTGACAGAATTCTCGACGCAGTCACAGCAAAAAACCGAAGGCGGCTAGAAATCAACAACTTTCCAGTTCGAGTGTTTCTAAATTCTGGTTATCCCGTTGAGCAAATATCCTAAATTTTAAGCAAAATGGCTGGACGCGATGCGGCTTCCAATCAGTTGATTGACTACAAAAACTCCCAAACGGTGAGTGTGTTGGTGTGCGAGCGCGAGAGAGCGCAACTGCATGTGTGTGCGTGTGTGTGCCGGCAGGAAAATAAGAAAACAAGAAAGGCAAAAAAGACAAAAAGAGAAGGGAGAAGCATTGTACGAAATGAAATAAACAAAAAACAAGACAATTTCAGTAACGATTGCTATGAGCCTGTACTTTCCATTAATTCCATTAATTATTAAAATGCGAAAAACAACAAAGGCTTGTCCAACAAATGCCAGTAATG$"
gene = "ACAAGTACGTGGTGCTGGGCTATAAAAACAAATGGAACTGCGCTGCAGTTGCTTTAGTTGACAGAATTCTCGACGCAGTCACAGCAAAAAACCGAAGGCGGCTAGAAATCAACAACTTTCCAGTTCGAGTGTTTCTAAATTCTGGTTATCCCGTTGAGCAAATATCCTAAATTTTAAGCAAAATGGCTGGACGCGATGCGGCTTCCAATCAGTTGATTGACTACAAAAACTCCCAAACGGTGAGTGTGTTGGTGTGCGAGCGCGAGAGAGCGCAACTGCATGTGTGTGCGTGTGTGTGCCGGCAGGAAAATAAGAAAACAAGAAAGGCAAAAAAGACAAAAAGAGAAGGGAGAAGCATTGTACGAAATGAAATAAACAAAAAACAAGACAATTTCAGTAACGATTGCTATGAGCCTGTACTTTCCATTAATTCCATTAATTATTAAAATGCGAAAAACAACAAAGGCTTGTCCAACAAATGCCAGTAATGATGAAGTCGTTTCCCCAGAGTTCGAAGTTCAAAAGCCGCCTCTCTGGGTCTTCTATTTATTTGTAGAAATCTTATCAGAAGTGAACTGATAAACCTTATCAGTGTGGCAATTCTTATTGGATTTGTGTTTTCTGTACCTAATCCAAAGTTAATTTTATTTTCCCTTTTATTTTTATTTTCCTTCTTAAGGAAATCCAAAGTTTAGTCTAGCTAACCTATTTTCTAACCTTGGGCGACCTCTTACCTTAGTCTATTTTCATTTTTCAGCGCTGTTTTCACTTAGCACTTCCAATTTCATTTGGCACTTGCATTGTGCATTGGATTTGAGTATTCGAATTTGTTCGTTTGCACTTTTGTGTTAATTATTCGCTTTCGCTTTTGGTTTAGATTAAATGCATGATTTCAATTGTGGAATCGTGCCCTCTTAAAATTTTGAACAAATATTGTTTGTTCACACTCGGAAAATTGTACTTTCCTAAATTCGTGGCATTCCACAATTGCAACATTTTATACACTGCAATCATTTTTGAAATCAAATATGCTGTTACATTTTATAAAAAGGAATGAAATATAACAAAATCAATAGTAAGCCTCACCCATAGATACTGTATCTATATATAATTTACTCGGCAATTTCCCTGTGGTTTTCCCAATGCCCTCTAAAAATCGATAAATGAAATAATAATTCTGCATAGCTATTTGTCAGCCAGTTGATAGGCATCCAAAACAGCTAAGTAATTGGCCCAAATAGATTAGATGTATCAGCCGTATAGATATTGTTGAGTGCCCACATTTGGTTAGTAACTAAAGTAAATTTCCATTGGATTTCGATACTTTGTAAGTATGAATTGTCAGTTCTCAAAGAATACAATCATTTTTTGTTTTTTAATTTAGTATACCTTTTTTTTACCCTATTCAATTAAGGCTAAATCTAGAACTTAGTTCTCGCATCAGCGATTTAAGTGCATTTGTGCATTGAGTCACGAAGCCAATTACTCCACTAACTAAAAGCAATCAGCGGTCATTATTTTATATAAACAAATATAGATATGAACTTTTTGAGAAACGCGATCGCTGGATTGCACCGCTTGTAAGTTGGAAAATTCGTTGGGAGTCAATGGACTCTGCGGCGACCGGCTCAATGTTGCACAAACTGGTTGCAATGCGCATCGACTGCGAGTGATTCGAGTGACAACGAACCCGGCTTTAAACCCGTTTTGGGTTGAGTGCGGAAGGCTAGACCGCTTCGGAAAATGCGTGTGCCAGTGCCATAAGTCAACGCAAATTTTCTTTCTCCAATTTGGTTATGCAACCTGGCTGGAACGCCCGGCTAAAAGCGACTGTATGTGAGAATATGTTTCGAACTGACCCAGAAAAGAGAGCTAGCTCCACAGTTGGTTGGTTGGTTCGTGTTGTTTTAGTCGATTTTGACCTTCTATGACGAAACTGTGAAAACTGCATTCTGTCAGGCAATTAAATGTTTGCCTTTGTTTCGGTGAATAATCATGGCAACTGGGTTCTTTGTTCTTTTGGGAATCTCGCAATCGAGTTAATTTGTTATCCAAGAAATTCGAAAATGCAAGGAAAATACTTGGAATATTTTGCAACGCAATGCATATTGCGTAAGGTCATTGTGAGTCGGTTAAGAACCCTATAAAAATCCTTATTTAAAATCCTTTAGCTTATTTTCTATGCATTCTCATTCTCGAATGAGTAACCTTCTGCTGATAAATTATGTCTAGAAGGCAAATTTATCTCATACTATCAGGTGCACTGATCTTATAAAGCCGCCTCATTAAGTTGCTGTTTGCCATCTTTGTGTCATTGGCATCGCCAACAAGTTGCAATTAAAAACTTTCCTGTAGACTGTGTTCTTTGCATTTCTATTCGAGATCGGTTGCGGTTCTCTTATCAGTTTGCTCTTCGCTGGAGATATCTACAACCAAATACATACATACATATATAAGAACGAAGGCTCTGATAGGAACCAGATTGACTTGTTGGACGTGGATCGCAGTTCCAGCGGTTCCACTCGCATTGTCATTCCCGGTTTGCTTCACTTATCGTATCTGCCCCCTGAAGTGGCCCCATATTAATGTCTCAAAGTGCAAAACATTCCGGCTGAATCGGCCGAATCTAAAGGGTGACTCCAGTATAAATTCCAGAATGATAATTTCGAAATTATATTCATTTCTTCGTTTTTATGTAGCGCTTTTTGTTTGTCAATAAGTGGGTATACTACTGGCCTGTCTAAGCTGATTTTATTATTAAATAGATACATATATACTTCTAGGAACTAAGATTTGATAAAGTCAACTACTAATTCCTTAAACAATTGTAAAAGTCTTAGTGCCCAACCGGGTAATCCCAACCCAAGTGCTATATGTATCTAAATGACTTTGTTGTGAAATATAAACAGTAGCTTAATTTATATGTTTGACACACCACCACCACCACCACCGGTGGCTTCACTTTGATAAACAACCTTCTGTATAAAATTCTATATGCATTGTTCTAGAAATGCTTTCGAATTAGAATTTCGCATGGCGTTGTTCAGGTTTCAAATACATTAGCTAAACAGTGCGCATTGGCCGGTGACAGTTTTGCAAATTGTTGGCTGAGTTAATTTTATAATGAGTGGCTGGCCTGCATGCGAAATAAACGATAAGGTGATTAAGGAGCGGTCTCGGATTACTCTAGGGAAATGAATTGAAAACAAGGAGTATGATTAGATCGACTGATAAATGACAGGACAGCTGGTGGGCTAATTGATATTTCAGATGGACAATACACTTTTATTAGCTTAATGGAATGAAATAATTAATTGATCGCCAGCTAATAAACTGTGTAACATCATTTTGATGTTTGAATACTTTCATTTAATTGGTTTTTTCCTTGCACACGGATTTGATATCTTGACATGTTAATGCAAAAAACCCATACATCATCACTATATTAAGGCCTCAAGCTCTTAAGAACACCTGGGATCTCTTATTTTAAATTTATATGCACTGCAAAATGACCCCCCGTGTAGAGATAACGCTCTGCAATCTTTATGTTAATGTTTTTTTCTTTTGGTGAAAATCGAATGATTAAATGAAAGTCCGCGGGCCAAAAATGAATTATCGCAATAGGCCCGATAAGATTGTAAGACACTTGAACCAATATATTGTTGCTGAATGTTGCATAGCATTTTAAAATTAAAATACATTTGATTAAACTTGACTTTCATTTGCTTTTGCAGGTTTCTCCTGGTGCCATTACCACCGGCAATGGAGCACCCATTGGAATCAAGGATGCCTCCCAGACGGTTGGTCCGCGAGGACCTATCCTGCTGCAGGATGTGAACTTCCTGGATGAGATGTCGCACTTCGACAGGGAGCGGATTCCAGAGCGTGTCGTGCACGCCAAGGGAGCTGGTGCTTTTGGTTACTTTGAGGTGACCCATGACATCACCCAGTATTGTGCCGCCAAGATTTTCGACAAGGTCAAGAAGCGCACTCCACTGGCCGTGCGATTCTCCACCGTGGGTGGTGAGAGCGGATCTGCGGACACCGCCCGCGATCCTCGAGGATTTGCCGTCAAGTTCTACACTGAGGATGGCGTCTGGGATTTGGTTGGCAACAACACGCCGGTCTTCTTCATTCGCGACCCGATCCTATTCCCCAGCTTCATTCACACCCAGAAGCGCAACCCGCAGACGCATCTGAAGGATCCGGACATGTTCTGGGACTTCCTCACCCTGCGACCGGAGTCCGCTCACCAGGTGTGCATCCTGTTCAGCGATCGCGGCACCCCCGACGGTTACTGCCACATGAACGGCTATGGCTCGCACACCTTCAAATTGATCAACGCCAAGGGCGAGCCCATCTATGCCAAGTTCCACTTCAAGACGGACCAGGGCATCAAGAATCTGGACGTGAAGACCGCCGATCAGTTGGCTAGCACTGATCCGGATTACAGCATTCGCGATCTGTACAACAGGATCAAGACCTGCAAGTTCCCCAGTTGGACGATGTACATTCAGGTCATGACCTACGAGCAGGCCAAGAAGTTCAAGTACAACCCCTTCGATGTCACCAAGGTCTGGTCGCAGAAGGAGTACCCTCTGATTCCTGTGGGCAAAATGGTGCTGGATCGCAATCCCAAGAACTACTTTGCTGAGGTTCGTTTTTTTTTTTTGTTTCAAATAGATTCTTCTTCATTTGAGATTAATGCAATCAAAAGAAAGATAACTAATATGGTAATTTTCTGACCCTCAGGTGGAGCAGATCGCCTTCAGTCCCGCTCACCTGGTGCCCGGCGTTGAGCCCTCTCCGGACAAGATGCTGCATGGTCGTCTGTTCTCCTACTCGGACACCCATCGCCATCGCCTGGGACCGAACTACTTGCAGATCCCGGTGAACTGCCCGTACAAGGTGAAGATTGAGAACTTCCAGCGGGATGGAGCCATGAATGTGACGGACAACCAGGATGGTGCCCCCAACTACTTCCCCAACTCGTTCAACGGTCCCCAGGAATGCCCCAGGGCCAGGGCCTTGTCCTCCTGCTGTCCGGTGACTGGAGATGTCTACCGCTACAGCAGCGGCGACACCGAGGACAACTTCGGCCAGGTCACCGACTTCTGGGTGCATGTGCTCGACAAGTGCGCCAAGAAGCGTCTGGTGCAGAACATTGCCGGCCATTTGAGCAACGCCAGCCAGTTCTTGCAGGAGCGGGCCGTCAAGAACTTCACCCAGGTGCACGCCGATTTCGGTCGCATGCTGACCGAGGAACTCAACCTGGCCAAGTCCTCGAAGTTCTAAGCTGAGCGAGCGGATTCGACGGATCAGACTTGGTTTTTGGCCTTAATTACGATTAATCCAATGGAACTAATTATTCCAACACCAACACCACCACCACCAACACCACCCATTCCGAAATTGACTAAGAAAGGCGAACGCTTTTCGGATCTTTGGTTGCGGGCTTAAATTGATTTAACCACGAAATGTGTCTTGTCTTGTTTTGATATTCAAAATTGTTTAAGTCATAAACAAATTACACTTACAATGTATATGTATATGCAATGAATATATAATAAATTCTTTATTTTTTG$"
gene2 = gene1[:len(gene1)-1] + gene
def sufar(txt):

    txt = txt + chr(0)
    N, tokens = len(txt), sorted(set(t for t in txt))

    equivalence = {t: i for i, t in enumerate(tokens)}
    N = len(txt)
    #equivalence = {'$':0, 'A':1, 'C':2, 'T':3, 'G':4}
    c, r = [equivalence[t] for t in txt], [(0, 0, 0)]

    for i in range(1, ceil(log2(N)) + 1):
        n = 2 ** (i - 1)

        r = sorted([(c[j],  c[(j + n) % N], j) for j in range(N)])

        c[r[0][2]] = 0
        for j in range(1, N):
            c[r[j][2]] = c[r[j - 1][2]]
            if r[j][0:2] != r[j - 1][0:2]:
                c[r[j][2]] += 1
    
    return [result[2] for result in r][1:]

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

    """def sortSuffixes(i, j, pos):

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

    return list(A)"""
    A = [len(s) - 1] + sufar(s[:len(s) - 1])
    return A
    
    pass

A = get_suffix_array(gene1)
print(A)


# In[81]:


def get_bwt(s, sa):
    """
    Input:
        s: a string terminated by a unique delimiter '$'
        sa: the suffix array of s

    Output:
        L: BWT of s as a string
    """
    L = ""
    for i in range(len(sa)):
        L += s[sa[i] - 1 % len(s)]
    return L

L = get_bwt(gene1, A)
print(L)


# In[82]:


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

print(get_F(L))


# In[83]:


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


# In[67]:


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

occ = get_occ(L)


# In[68]:


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
        prev_sp, prev_ep, = sp, ep
        pos -= 1
        c = p[pos]
        sp = M[c] + occ[c][sp - 1]
        ep = M[c] + occ[c][ep] - 1
        if sp > ep:
            return ((prev_sp, prev_ep), len(p) - pos - 1)
    
    
    return ((sp, ep), len(p))
        


# In[69]:


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
        while len(rest) > 0:
            result = exact_suffix_matches(rest, M, occ)
            print("string:", rest)
            print("max suffix match:", rest[len(rest) - result[1]:])
            print(result)
            print("--------------")
            start = result[0][0]
            end = result[0][1]
            if end - start <= 5:
                for i in range(start, end + 1):
                    hit = (len(rest) - result[1], A[i], result[1])
                    #(position in read, position in genome, length)
                    hits.append(hit)
            rest = rest[:len(rest) - result[1]]
            
        hits.sort(key=lambda x: x[1])
        print(hits)
        i = 0
        while i < len(hits):
            print("considering hit:", hits[i])
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
                print("extending:", hits[curr])
                genomic_pos = hits[curr][1] + hits[curr][2]
                for j in range(1, hits[next_hit][1] - hits[curr][1]):
                    if r[end + j] != T[genomic_pos + j]:
                        mismatches += 1
                end = hits[next_hit][0] + hits[next_hit][2]
                curr = next_hit
                next_hit += 1
                num_extended += 1
                print(next_hit)
            i += num_extended
            print(mismatches)
            while end < len(r) - 1:
                end += 1
                if r[end] != T[genomic_start + end]:
                    mismatches += 1    
            print(mismatches)
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
        
                        
            


# In[77]:


p =  "TCAGTTGATTGACTACAAAA"
p1 = "CAAAAAGAGAAGGGAGAAGCA"
p2 = "CGTGGTGCTGGGCTATAAA"
p3 = "AAATGCCAGTAATG"
p4 = "AAATGGCTGGACGCGATGCGG"
q =  "TCAGTGGATTGATTCCAAAA"
q1 = "CATAAAGAGAAGGGAGGGGCA"
q2 = "CGTGGTCCTGGGCTATAGA"
q3 = "ATATGCCAGTAATA"
q4 = "AAATGGCTGGACGCGACGCGG"


A = get_suffix_array(gene)
L = get_bwt(gene, A)
F = get_F(L)
M = get_M(F)
occ = get_occ(L)


# In[78]:


gdi = exact_suffix_matches(p, M, occ)
print(A[gdi[0][0]])
print(gdi)


# In[79]:


reads = [q, q1, q2, q3, q4]
alignments = simple_align(reads, gene)


# In[73]:


print(alignments)


# In[245]:


p2 = "TAAGTTGATTGACT"
gdi2 = exact_suffix_matches(p2, M, occ)
print(gdi2)
print(A[gdi2[0][0]])


# In[145]:

def find_sp(lo, hi, c, pos, SA, s):

    if(lo == hi):
        if (s[SA[lo]+pos] == c):
            return lo
        else:
            return None

    while lo < hi - 1:

        mid = (lo + hi)//2
        if (s[SA[mid] + pos] < c):
            lo = mid
        else:
            hi = mid

    if (lo == hi - 1):
        if (s[SA[lo] + pos] == c):
            return lo
        elif (s[SA[hi] + pos] == c):
            return hi
        else:
            return None

    return lo + 1

def find_ep(lo, hi, c, pos, SA, s):

    if(lo == hi):
        if (s[SA[lo]+pos] == c):
            return lo
        else:
            return None

    while lo < hi - 1:

        mid = (lo + hi)//2
        if (s[SA[mid] + pos] > c):
            hi = mid
        else:
            lo = mid

    if (lo == hi - 1):
        if (s[SA[lo] + pos] == c) and (s[SA[hi] + pos] != c):
            return lo
        elif (s[SA[hi] + pos] == c):
            return hi
        else:
            return None

    return hi - 1


def binary_search(seq, SA, s):

    sp = 0
    ep = len(s) - 1
    pos = 0

    while sp <= ep and pos < len(seq):

        prev_sp, prev_ep = sp, ep
        c = seq[pos]

        if (sp == ep):
            if (s[SA[sp] + pos] == c):
                pos = pos + 1
                continue
            else:
                return (sp, ep), pos

        sp = find_sp(prev_sp, prev_ep, c, pos, SA, s)
        ep = find_ep(prev_sp, prev_ep, c, pos, SA, s)
        if (sp == None) or (ep == None):
            return (prev_sp, prev_ep), pos
        pos = pos + 1

    if (sp > ep):
        return (prev_sp, prev_ep), pos
    elif (sp == None) or (ep == None):
        return (prev_sp, prev_ep), pos-1
    else:
        return (sp, ep), pos



def simple_sa_align(reads, T):
    A = get_suffix_array(T)
    alignments = []
    print("reee?")
    for r in reads:
        rest = r
        potential_aligns = []
        hits = []
        while len(rest) > 0:
            result = binary_search(rest, A, T)
            print("string:", rest)
            print("max prefix match:", rest[:result[1]])
            print(result)
            print("--------------")
            start = result[0][0]
            end = result[0][1]
            if end - start <= 5:
                for i in range(start, end + 1):
                    hit = (len(r) - len(rest), A[i], result[1])
                    #(position in read, position in genome, length)
                    print(hit)
                    hits.append(hit)
            rest = rest[result[1]:]
        lowest_mismatch = 7
        best = ()
        found_best = False
        print("hits ", hits)
        for h in hits:
            mismatches = 0
            start = h[0]
            genomic_start = h[1]
            end = h[0] + h[2] - 1
            while start > 0:
                start -= 1
                genomic_start -= 1
                if r[start] != T[genomic_start]:
                    mismatches += 1
                if mismatches > 6:
                    break
            if mismatches > 6:
                continue
            while end < len(r) - 1:
                end += 1
                if r[end] != T[genomic_start + end]:
                    mismatches += 1
                if mismatches > 6:
                    break
            if mismatches > 6:
                continue

            if mismatches == 0:
                alignments.append((start, genomic_start, len(r)))
                found_best = True
                break
            if mismatches < lowest_mismatch:
                best = (start, genomic_start, len(r))
                lowest_mismatch = mismatches
        if not found_best and best:
            alignments.append(best)
    return alignments

        
        
        





def sufar(txt):

    txt = txt + chr(0)
    N, tokens = len(txt), sorted(set(t for t in txt))

    equivalence = {t: i for i, t in enumerate(tokens)}
    N = len(txt)
    #equivalence = {'$':0, 'A':1, 'C':2, 'T':3, 'G':4}
    c, r = [equivalence[t] for t in txt], [(0, 0, 0)]

    for i in range(1, ceil(log2(N)) + 1):
        n = 2 ** (i - 1)

        r = sorted([(c[j],  c[(j + n) % N], j) for j in range(N)])

        c[r[0][2]] = 0
        for j in range(1, N):
            c[r[j][2]] = c[r[j - 1][2]]
            if r[j][0:2] != r[j - 1][0:2]:
                c[r[j][2]] += 1
    
    return [result[2] for result in r][1:]


#


def align_genome(reads, T, intron_lb, intron_ub):
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
    
    def extend_alignment_to_end(partial_alignment, read):
        potential_alignment = partial_alignment[0]
        mismatches = partial_alignment[1]
        hit = potential_alignment[len(potential_alignment) - 1]
        new_len = hit[2]
        while hit[0] + new_len < len(read):
            if read[hit[0] + new_len] != T[hit[1] + new_len]:
                mismatches += 1   
            if mismatches > 6:
                return None
            new_len += 1
        new_hit = (hit[0], hit[1], new_len)
        return (potential_alignment[:len(potential_alignment) - 1] + [new_hit], mismatches)

        
    
    def extend_alignment(partial_alignment, last_hit, read):
        potential_alignment = partial_alignment[0]
        print(partial_alignment)
        mismatches = partial_alignment[1]
        first_hit = potential_alignment[len(potential_alignment) - 1]
        dist = last_hit[0] - first_hit[0] - first_hit[2]
        print(first_hit, last_hit)
        intron_size = last_hit[1] - first_hit[1] - first_hit[2] - dist
        best_cost = 0
        best_i = 0
        for i in range(1, dist + 1):
            if read[last_hit[0] - i] != T[last_hit[1] - i]:
                best_cost += 1
        M = np.zeros(dist + 1, dtype=int)
        M[0] = best_cost
        for i in range(1, dist + 1):
            M[i] = M[i-1]
            #print(read[first_hit[0] + first_hit[2] - 1 + i] + " " + T[first_hit[1] + first_hit[2] - 1 + i])
            #print(read[first_hit[0] + first_hit[2] - 1 + i] + " " + T[first_hit[1] + first_hit[2] - 1 + intron_size + i])
            #print("-----------------------")
            if read[first_hit[0] + first_hit[2] - 1 + i] != T[first_hit[1] + first_hit[2] - 1 + i]:
                M[i] += 1
            if read[first_hit[0] + first_hit[2] -1 + i] != T[first_hit[1] + first_hit[2] - 1 + intron_size + i]:
                M[i] -= 1
            if M[i] < best_cost:
                best_cost = M[i]
                best_i = i
        print(M)
        if mismatches + best_cost > 6:
            return None
        new_left = (first_hit[0], first_hit[1], first_hit[2] + best_i)
        dl = dist - best_i
        new_right = (last_hit[0] - dl, last_hit[1] - dl, last_hit[2] + dl)
        return (potential_alignment[:len(potential_alignment) - 1] + [new_left, new_right], mismatches + best_cost)
    
    
    def extend_right(partial_alignment, i, hits, read):
        potential_alignment = partial_alignment[0]
        mismatches = partial_alignment[1]
        valid_aligns = []
        a = extend_alignment_to_end(partial_alignment, read)
        if a:
            valid_aligns.append(a)
        if len(potential_alignment) < 3:
            print(potential_alignment)
            print(i)
            for j in range(i + 1, len(hits)):
                hit = hits[i]
                next_hit = hits[j]
                genomic_dist = next_hit[1] - hit[1] - hit[2]
                if genomic_dist > intron_ub:
                    break
                dist = next_hit[0] - hit[0] - hit[2]
                if dist < 0:
                    continue
                intron_length = genomic_dist - dist
                if intron_length >= intron_lb and intron_length <= intron_ub:
                    pa = extend_alignment(partial_alignment, next_hit, read)
                    print("pa is " + str(pa))
                    print(j)
                    if pa:
                        valid_aligns.extend(extend_right(pa, j, hits, read))
        return valid_aligns

    for r in reads:
        rest = r
        complete_aligns = []
        hits = []
        while len(rest) > 0:
            result = exact_suffix_matches(rest, M, occ)
            #print("string:", rest)
            #print("max suffix match:", rest[len(rest) - result[1]:])
            #print(result)
            #print("--------------")
            start = result[0][0]
            end = result[0][1]
            if end - start <= 5:
                for i in range(start, end + 1):
                    hit = (len(rest) - result[1], A[i], result[1])
                    #(position in read, position in genome, length)
                    if len(r) - hit[0] <= len(T) - hit[1]:
                        hits.append(hit)
            rest = rest[:len(rest) - result[1]]
            
        hits.sort(key=lambda x: x[1])
        print(hits)
        for i in range(len(hits)):
            start = hits[i][0]
            genomic_start = hits[i][1]
            mismatches = 0
            len_so_far = hits[i][2]
            #extend left
            while start > 0:
                start -= 1
                genomic_start -= 1
                len_so_far += 1
                if r[start] != T[genomic_start]:
                    mismatches += 1
                if mismatches > 6:
                    break
            if mismatches > 6:
                i += 1
                continue
            partial_alignment = ([(0, genomic_start, hits[i][2] + hits[i][0] - start)], mismatches)
            v = extend_right(partial_alignment, i, hits, r)
            complete_aligns.extend(v)
            
        best = ()
        min_mismatches = 7
        for x in complete_aligns:
            if x[1] < min_mismatches:
                best = x
        alignments.append(best[0])
        
    return alignments


# In[375]:



if __name__ == "__main__":
    """p = "GCGGCTTCCAATCAGTTGATTGCGTTTCCCCAGAGTTC"
    p1 ="GCGGCTTGCAATCAGTTGATTGCGTTTCCCCAGAGTTC"
    p2 ="GCAACTACCAATCAGTAGATTACGATTCCCGACAGTTC"""
    A = get_suffix_array(gene)
    p = "TAAACAAAAAACAAGACAATTTCAGTAACGATTGCTATGAGCCT"
    p2 ="TAGACAAGAATCAAGACAATTTCAGTGACAATTGCTGTGAGCCT"
    p3 ="TAAACAAGAAACAAGACAATTTCAGTTTATTGGATTTGT"
    m = binary_search(p, A, gene)
    print(m)
    print("genomic pos ", A[m[0][0]])
    reads = [p3]
    a = align_genome(reads, gene, 20, 500)
    print(a)





