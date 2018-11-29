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

from project import *
from math import ceil, log2
import numpy as np
import time


ALPHABET = [TERMINATOR] + BASES


def to_int_keys_best(l):
    """
    l: iterable of keys
    returns: a list with integer keys
    """
    seen = set()
    ls = []
    for e in l:
        if not e in seen:
            ls.append(e)
            seen.add(e)
    ls.sort()
    index = {v: i for i, v in enumerate(ls)}
    return [index[v] for v in l]

def get_suffix_array(s):

    if not s:
        return []
    s += chr(0)

    equivalence = {t: i for i, t in enumerate(sorted(set(s)))}
    cls = np.array([equivalence[t] for t in s])
    ns = 2**np.arange(ceil(log2(len(s))))

    for n in ns[:-1]:
        cls1 = np.roll(cls, -n)
        inds = np.lexsort((cls1, cls))
        result = np.logical_or(np.diff(cls[inds]), 
                               np.diff(cls1[inds]))

        cls[inds[0]] = 0
        cls[inds[1:]] = np.cumsum(result)

<<<<<<< HEAD
    cls1 = np.roll(cls, ns[-1])

    return np.lexsort((cls1, cls))[1:]
=======
    return list(A)
>>>>>>> 32fb04f9a4604ddb28ca8f45a0f3fa4add685019


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
        prev_sp, prev_ep, = sp, ep
        pos -= 1
        c = p[pos]
        sp = M[c] + occ[c][sp - 1]
        ep = M[c] + occ[c][ep] - 1
        if (ep >= len(occ)) :
            ep -= 1
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
        genome_loc = {}
        trans_ix = 0
        genome_loc_array = []

        for gene in known_genes:
            for isoform in gene.isoforms:
                iso = ""

                for exon in isoform.exons:
                    iso = iso + genome_sequence[exon.start:exon.end]
                    exon_len = len(genome_sequence[exon.start:exon.end])
                    keeper = [exon.start, len(genome_sequence[exon.start:exon.end])]
                    genome_loc[trans_ix] = tuple(keeper)
                    genome_loc_array[trans_ix:trans_ix+exon_len] = list(range(exon.start,exon.end))
                    trans_ix = trans_ix + len(iso)

                transcriptome = transcriptome + iso + "$$$$$$"
                genome_loc_array[trans_ix:trans_ix+6] = [-1]*6
                trans_ix = trans_ix + 6

        self.genome_loc = genome_loc
        self.trans_SA = get_suffix_array(transcriptome)
        self.transcriptome = transcriptome
        self.genome_loc_array = genome_loc_array
        self.genome_sequence = genome_sequence
        self.genome_SA = []
        self.genome_SA = get_suffix_array(self.genome_sequence)
        self.L = get_bwt(self.genome_sequence, self.genome_SA)
        self.F = get_F(self.L)
        self.M = get_M(self.F)
        self.occ = get_occ(self.L)
        self.intron_lb = 20
        self.intron_ub = 10000

        pass


    def index_translator(self, aligned_piece):
    
        def trans_to_genome(i):

            ix = self.genome_loc
            offset = 0

            while i >= 0:
                if i in ix:
                    return ix.get(i)[0] + offset
                else:
                    i = i - 1
                    offset = offset + 1

        read_start = aligned_piece[0]
        genome_start = trans_to_genome(aligned_piece[1])
        read_len = aligned_piece[2]

        return (read_start, genome_start, read_len)

    def find_sp(self, lo, hi, c, pos, SA, s):

        if (lo == hi):
            if (s[SA[lo] + pos] == c):
                return lo
            else:
                return None

        while lo < hi - 1:

            mid = (lo + hi) // 2
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

    def find_ep(self, lo, hi, c, pos, SA, s):

        if (lo == hi):
            if (s[SA[lo] + pos] == c):
                return lo
            else:
                return None

        while lo < hi - 1:

            mid = (lo + hi) // 2
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

    def binary_search(self, seq, SA, s):

        sp = 0
        ep = len(s) - 1
        pos = 0

        while sp <= ep and pos < len(seq):

            prev_sp = sp
            prev_ep = ep
            c = seq[pos]

            if (sp == ep):
                if (s[SA[sp] + pos] == c):
                    pos = pos + 1
                    continue
                else:
                    return (sp, ep), pos

            sp = self.find_sp(prev_sp, prev_ep, c, pos, SA, s)
            ep = self.find_ep(prev_sp, prev_ep, c, pos, SA, s)

            if (sp == None) or (ep == None):
                return (prev_sp, prev_ep), pos
            pos = pos + 1

        if (sp > ep):
            return (prev_sp, prev_ep), pos
        elif (sp == None) or (ep == None):
            return (prev_sp, prev_ep), pos-1
        else:
            return (sp, ep), pos
        
    def extend_alignment_to_end(self, partial_alignment, read):
        T = self.genome_sequence
        intron_lb = self.intron_lb
        intron_ub = self.intron_ub
        A = self.genome_SA
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
    
    
    def extend_alignment(self, partial_alignment, last_hit, read):
        T = self.genome_sequence
        intron_lb = self.intron_lb
        intron_ub = self.intron_ub
        A = self.genome_SA
        potential_alignment = partial_alignment[0]
        mismatches = partial_alignment[1]
        first_hit = potential_alignment[len(potential_alignment) - 1]
        dist = last_hit[0] - first_hit[0] - first_hit[2]
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

        if mismatches + best_cost > 6:
            return None
        new_left = (first_hit[0], first_hit[1], first_hit[2] + best_i)
        dl = dist - best_i
        new_right = (last_hit[0] - dl, last_hit[1] - dl, last_hit[2] + dl)
        return (potential_alignment[:len(potential_alignment) - 1] + [new_left, new_right], mismatches + best_cost)
    

    def extend_right(self, partial_alignment, i, hits, read):
        T = self.genome_sequence
        intron_lb = self.intron_lb
        intron_ub = self.intron_ub
        A = self.genome_SA
        potential_alignment = partial_alignment[0]
        mismatches = partial_alignment[1]
        valid_aligns = []
        a = self.extend_alignment_to_end(partial_alignment, read)
        if a:
            valid_aligns.append(a)
        if len(potential_alignment) < 3:

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
                    pa = self.extend_alignment(partial_alignment, next_hit, read)

                    if pa:
                        valid_aligns.extend(self.extend_right(pa, j, hits, read))
        return valid_aligns


    def align_genome(reads, T, intron_lb, intron_ub, A):
        #take longest hits, try to extend to hits left and right
        #merge with nearby hits
        #once merged as much as possible, calculate score
        #Take highest scoring match
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
            mismatches = partial_alignment[1]
            first_hit = potential_alignment[len(potential_alignment) - 1]
            dist = last_hit[0] - first_hit[0] - first_hit[2]
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
                if read[first_hit[0] + first_hit[2] - 1 + i] != T[first_hit[1] + first_hit[2] - 1 + i]:
                    M[i] += 1
                if read[first_hit[0] + first_hit[2] -1 + i] != T[first_hit[1] + first_hit[2] - 1 + intron_size + i]:
                    M[i] -= 1
                if M[i] < best_cost:
                    best_cost = M[i]
                    best_i = i
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
            # print(hits)
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


    def index_translator(self, aligned_pieces):

        read_start = aligned_pieces[0]
        ref_start = aligned_pieces[1]
        len_match = aligned_pieces[2]
        ix = self.genome_loc_array
        alignments = []

        k = 1
        genomic_start = ix[ref_start]

        for i in range(ref_start + 1, ref_start + len_match):
            pos = ix[i]
            if ix[i] != ix[i-1] + 1:
                alignments.append((read_start, genomic_start, k))
                read_start += k
                genomic_start = pos
                k = 0
            k += 1

        alignments.append((read_start, genomic_start, k))

        return alignments
    
    def align_genome_single(self, r):

        T = self.genome_sequence
        intron_lb = self.intron_lb
        intron_ub = self.intron_ub
        A = self.genome_SA
        rest = r
        complete_aligns = []
        hits = []
        while len(rest) > 0:
            result = exact_suffix_matches(rest, self.M, self.occ)
            #print("string:", rest)
            #print("max suffix match:", rest[len(rest) - result[1]:])
            #print(result)
            #print("--------------")
            start = result[0][0]
            end = result[0][1]
            if end - start <= 8:
                for i in range(start, end + 1):
                    hit = (len(rest) - result[1], A[i], result[1])
                    #(position in read, position in genome, length)
                    if len(r) - hit[0] <= len(T) - hit[1]:
                        hits.append(hit)
            rest = rest[:len(rest) - result[1]]

        hits.sort(key=lambda x: x[1])

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
            v = self.extend_right(partial_alignment, i, hits, r)
            complete_aligns.extend(v)

        best = ([], 7)
        min_mismatches = 7
        for x in complete_aligns:
            if x[1] < min_mismatches:
                best = x
        return best[0]

    
    def simple_sa_align_single(self, r):
        A = self.trans_SA
        T = self.transcriptome
        rest = r
        potential_aligns = []
        hits = []
        while len(rest) > 0:
            result = self.binary_search(rest, A, T)
            """print("string:", rest)
            print("max prefix match:", rest[:result[1]])
            print(result)
            print("--------------")"""
            start = result[0][0]
            end = result[0][1]
            if end - start <= 5:
                for i in range(start, end + 1):
                    hit = (len(r) - len(rest), A[i], result[1])
                    #(position in read, position in genome, length)
                    hits.append(hit)
            rest = rest[result[1]:]
        lowest_mismatch = 7
        best = ()
        found_best = False
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
                return(start, genomic_start, len(r))
                found_best = True
                break
            if mismatches < lowest_mismatch:
                best = (start, genomic_start, len(r))
                lowest_mismatch = mismatches
        if not found_best and best:
            return best
        return []

<<<<<<< HEAD
=======

        pass
>>>>>>> 32fb04f9a4604ddb28ca8f45a0f3fa4add685019


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

        a1 = self.simple_sa_align_single(read_sequence)
        if a1:
            return self.index_translator(a1)
        
        a2 = self.align_genome_single(read_sequence)
        
        return a2
    
        pass
