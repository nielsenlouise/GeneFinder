# -*- coding: utf-8 -*-
"""
This code takes in a DNA sequence and returns the amino acid sequences that are most likely to be proteins.

@author: Louise Nielsen

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

        I added the second two tests to see if it works for everything.
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else: # in case one of the entries is not a nucleotide
        print "Input was not a nucleotide."



def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
        i tested an empty string. it works.
    >>> get_reverse_complement('')
    ''
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    comp_dna_list = map(get_complement,dna) # creates a list of the complementary strand of dna
    comp_dna = ''.join(comp_dna_list) # joins comp_dna_list into one string
    reverse_comp_dna = comp_dna[::-1] # reverses the string
    return reverse_comp_dna


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
        also test that it works if it doesn't end in a stop codon.
    >>> rest_of_ORF("ATGGATCTG")
    'ATGGATCTG'
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # stop codons: TAG, TAA, TGA
    for i in range(len(dna)/3):
        if dna[3*i:3*i+3] == 'TAG' or dna[3*i:3*i+3] == 'TAA' or dna[3*i:3*i+3] == 'TGA':
            return dna[:3*i]
    return dna

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
        it should only do things if things start with start codon
        it should also not let nested things happen
    >>> find_all_ORFs_oneframe("ATGATGTAGATCGATCGATCGATCGTAG")
    ['ATGATG']
    >>> find_all_ORFs_oneframe("CATGAATGTAGATAGATGTGCCC")
    ['ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    All_ORFs = []
    position = 0
    while position < len(dna):
        if dna[position:position + 3] == 'ATG':
            end_dna = len(rest_of_ORF(dna[position:]))
            All_ORFs.append(rest_of_ORF(dna[position:]))
            position += end_dna
        else:
            position += 3
    return All_ORFs




def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
        it returns an empty list when nothing happens
    >>> find_all_ORFs("TGAGCGTGCGTGCGTG")
    []
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    All_ORFs_all_frames = []
    All_ORFs_all_frames.extend(find_all_ORFs_oneframe(dna))
    All_ORFs_all_frames.extend(find_all_ORFs_oneframe(dna[1:]))
    All_ORFs_all_frames.extend(find_all_ORFs_oneframe(dna[2:]))
    return All_ORFs_all_frames



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
        what happens if there's nothing there
    >>> find_all_ORFs_both_strands("")
    []
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    All_ORFs_both_strands = []
    All_ORFs_both_strands.extend(find_all_ORFs(dna))
    All_ORFs_both_strands.extend(find_all_ORFs(get_reverse_complement(dna)))
    return All_ORFs_both_strands


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("")
    'longest_ORF received an empty string.'
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    if dna != '':
        return max(find_all_ORFs_both_strands(dna)) # finds and returns longest entry in list
    else:
        return 'longest_ORF received an empty string.'



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    Length_ORFs_noncoding = [] # set a list for the length of non-coding ORFs
    for i in range(num_trials): # for loop that runs for the number of random trials
        shuffled = shuffle_string(dna) # shuffles the input
        length_of_long_ORF = len(longest_ORF(shuffled)) # stores length of longest ORF
        Length_ORFs_noncoding.append(length_of_long_ORF) # stores lengths of ORFs in a list
#    print Length_ORFs_noncoding # prints the list of long ORF lengths
    return max(Length_ORFs_noncoding) # returns the longest of the lengths of ORFs


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    codon = ''
    sequence = ''
    for x in range(len(dna)/3):
        codon = dna[3*x:3*x+3]
        if len(codon) == 3:
            amino_acid = aa_table[codon]
        sequence += amino_acid
        x += 3
    return sequence



def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna,1500)
    All_ORFs = find_all_ORFs_both_strands(dna)
    Long_enough_ORFs = [i for i in All_ORFs if len(i) > threshold]
    Amino_acids = [coding_strand_to_AA(j) for j in Long_enough_ORFs]
    return Amino_acids




if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    #doctest.run_docstring_examples(coding_strand_to_AA, globals())
