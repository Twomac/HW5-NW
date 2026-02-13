# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    nw.align(seq1, seq2)

    # I worked out correct matrices by hand
    correct_align_matrix = np.array([[0, -10, -11, -12],
                                     [-10, 5, -5, -6],
                                     [-11, -5, 4, -6],
                                     [-12, -6, 0, 5],
                                     [-13, -7, -5, 5]])
    
    correct_gapA_matrix = np.array([[-np.inf, -10, -11, -12],
                                     [-np.inf, -20, -5, -6],
                                     [-np.inf, -21, -15, -6],
                                     [-np.inf, -22, -16, -10],
                                     [-np.inf, -23, -17, -15]])
    
    correct_gapB_matrix = np.array([[-np.inf, -np.inf, -np.inf, -np.inf],
                                     [-10, -20, -21, -22],
                                     [-11, -5, -15, -16],
                                     [-12, -6, -6, -16],
                                     [-13, -7, -7, -5]])
    
    assert np.array_equal(nw.align_matrix(), correct_align_matrix)
    assert np.array_equal(nw.gapA_matrix(), correct_gapA_matrix)
    assert np.array_equal(nw.gapB_matrix(), correct_gapB_matrix)

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    score, seqA_align, seqB_align = nw.align(seq3, seq4)
    
    # correct backtrace sequence is rather trivial, so I wrote it out by hand
    correct_seqA_align = "MAVHQLIRRP"
    correct_seqB_align = "M---QLIRHP"
    correct_score = 18

    assert seqA_align == correct_seqA_align
    assert seqB_align == correct_seqB_align
    assert score == correct_score




