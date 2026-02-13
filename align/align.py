# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        self._align_matrix = np.zeros((len(seqA) + 1, len(seqB) + 1)) # matrix for storing alignment scores
        self._gapA_matrix = np.zeros((len(seqA) + 1, len(seqB) + 1)) # matrix for storing scores of alignments that end with gap in A
        self._gapB_matrix = np.zeros((len(seqA) + 1, len(seqB) + 1)) # matrix for storing optimal scores of alignments that end with gap in B

        # Init matrices for backtrace procedure
        self._back = np.zeros((len(seqA) + 1, len(seqB) + 1))
        self._back_A = np.zeros((len(seqA) + 1, len(seqB) + 1))
        self._back_B = np.zeros((len(seqA) + 1, len(seqB) + 1))

        
        # TODO: Implement global alignment here
        '''
        Here we implement the Gotoh (1982) algorithm which is actually a 3 matrix approach to optimizing the global 
        alignment between two sequences. It's helpful for affine gap penalties because it allows us to track 
        whether we are opening or extending gaps. _align_matrix is still the matrix of scores for optimal 
        alignments M[i,j] with i residues from seqA and j residues from seqB. _gapA_matrix on the other hand
        stores the optimal scores for alignments between i and j (Q[i,j]) that END with a gap in seqA. _gapB_matcrix are
        optimal scores for alignments that end with gaps in seqB (P[i,j]). Since we have 3 matrices of scores, we also
        have 3 backtrace matrices because we will be hopping between all 3 matrices during our traceback.
        When we are in _align_matrix (for instance at the beginning of the traceback when we start at M[len(a),len(B)]),
        _back will tell us whether the alignment that we came from was the upper left diagonal in _align_matrix 
        (sequence that ends in no gaps) with a value of 0, or to move to the sequence at the same indices in _gapA_matrix
        (sequence that ends with gap in A at Q[i,j]) with a value of 1, or to move to the sequence at the same indices in
        _gapB_matrix (sequence that ends with gap in B at P[i,j]). _back_A will take a value of 0 if we should move to 
        the left column but back to _align_matrix and a value of 1 if we should stay in _gapA_matrix. Similar for _back_B
        but moving to the upper row instead of left column.
        '''
        # we begin with initializing the first row and column

        # upper left cell
        self._align_matrix[0,0] = 0 
        self._gapA_matrix[0,0] = -np.inf
        self._gapB_matrix[0,0] = -np.inf
        
        # top row, pairing residues in seqB with gaps in A
        for j in range(1, len(seqB) + 1):
            self._align_matrix[0,j] = self.gap_open + (j-1)*self.gap_extend  # opening and extending gaps in A for each residue in seqB
            self._gapA_matrix[0,j] = self.gap_open + (j-1)*self.gap_extend  # opening and extending gaps in A for each residue in seqB
            self._gapB_matrix[0,j] = -np.inf  # impossible to end in a gap in B if we have not iterated over any residue from seqA
            self._back[0,j] = 1 # we should move to _gapA_matrix to find alignments ending with gaps in A
            self._back_A[0,j] = 1 # stay in _gapA_matrix since we have consumed all residues in seqA
            self._back_B[0,j] = -1 # huh? we can't end in a gap in B if we are out of residues in seqA!

        # left column, pairing residues in seqA with gaps in B
        for i in range(1, len(seqA) + 1):
            self._align_matrix[i,0] = self.gap_open + (i-1)*self.gap_extend
            self._gapA_matrix[i,0] = float('-inf')
            self._gapB_matrix[i,0] = self.gap_open + (i-1)*self.gap_extend
            self._back[i,0] = 2 # we should move to _gapB_matrix to find alignments ending with gaps in B
            self._back_A[i,0] = -1 # huh? we can't end in a gap in A if we are out of residues in seqB!
            self._back_B[i,0] = 1 # stay in _gapB_matrix since we have consumed all residues in seqB

        # fill out rest of matrices
        for i in range(1, len(seqA) + 1):
            for j in range(1, len(seqB) + 1):
                
                # calculate optimal scoree for alignments that end in gaps in B (P[i,j])
                insert_Bgap = self._align_matrix[i-1,j] + self.gap_open # score for opening a gap in B paired with next residue in seqA (we look at score in cell above current indices in _align_matrix and add gap open penalty since we are only moving down one row and not consuming any residues in seqB)
                extend_Bgap = self._gapB_matrix[i-1,j] + self.gap_extend # score for extending a gap in B paired with next residue in seqA (since we are coming from the matrix of alignments that end with gaps in B, we know we are merely extending)
                self._gapB_matrix[i,j] = max(insert_Bgap, extend_Bgap) # optimal score for alignments that end with gap in B
                if self._gapB_matrix[i,j] == insert_Bgap:
                    self._back_B[i,j] = 0 # move to _align_matrix since that is where we came from when we opened the gap in B
                else:
                    self._back_B[i,j] = 1 # stay in _gapB_matrix since we extended a gap in B
                
                # calculate optimal score for alignments that end in gaps in A (Q[i,j])
                insert_Agap = self._align_matrix[i,j-1] + self.gap_open # score for opening a gap in A paired with next residue in seqB (we look at score in cell to the left of current indices in _align_matrix and add gap open penalty since we are only moving right one column and not consuming any residues in seqA)
                extend_Agap = self._gapA_matrix[i,j-1] + self.gap_extend # score for extending a gap in A paired with next residue in seqB (since we are coming from the matrix of alignments that end with gaps in A, we know we are merely extending)
                self._gapA_matrix[i,j] = max(insert_Agap, extend_Agap) # optimal score for alignments that end with gap in A
                if self._gapA_matrix[i,j] == insert_Agap:
                    self._back_A[i,j] = 0 # move to _align_matrix since that is where we came from when we opened the gap in A
                else:
                    self._back_A[i,j] = 1 # stay in _gapA_matrix since we extended a gap in A

                # calculate optimal score for alignments that end in no gaps (M[i,j])
                match = self._align_matrix[i-1,j-1] + self.sub_dict[(seqA[i-1], seqB[j-1])] # score for pairing next residues in both seqA and seqB
                self._align_matrix[i,j] = max(match, self._gapA_matrix[i,j], self._gapB_matrix[i,j]) # optimal score of all possible alignments at this cell (end in Agap, end in Bgap, or end in no gaps)
                if self._align_matrix[i,j] == match:
                    self._back[i,j] = 0 # move to upper left diagonal in _align_matrix since we came from an alignment that ended in no gaps
                elif self._align_matrix[i,j] == self._gapA_matrix[i,j]:
                    self._back[i,j] = 1 # move to _gapA_matrix since we came from an alignment that ended in a gap in A
                else:
                    self._back[i,j] = 2 # move to _gapB_matrix since we came from an alignment that ended in a gap in B	
        		    
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        positions = ['M', 'Q', 'P'] # keeping track of which matrix we are in during backtrace
        i = len(self._seqA) # start at bottom right cell of matrices
        j = len(self._seqB)
        ε = 'M' # start in _align_matrix since that is where we have the optimal score for the global alignment of the two sequences
        self.alignment_score = self._align_matrix[i,j] # optimal score for global alignment is in bottom right cell of _align_matrix
        reverse_A_align = '' # we will build the alignment in reverse order and then reverse it at the end
        reverse_B_align = ''

        while i > 0 or j > 0: # continue until we have reached the upper left cell of the matrices
            if ε == 'M': # we are in _align_matrix
                if self._back[i,j] == 0: # move to upper left diagonal in _align_matrix since we came from an alignment that ended in no gaps
                    reverse_A_align += self._seqA[i-1]
                    reverse_B_align += self._seqB[j-1]
                    i -= 1
                    j -= 1
                elif self._back[i,j] == 1: # move to _gapA_matrix since we came from an alignment that ended in a gap in A
                    ε = 'Q'
                else: # move to _gapB_matrix since we came from an alignment that ended in a gap in B	
                    ε = 'P'
            elif ε == 'Q': # we are in _gapA_matrix
                if self._back_A[i,j] == 0: # move to _align_matrix since that is where we came from when we opened the gap in A
                    ε = 'M'  # move back to _align_matrix
                reverse_A_align += '_' # gap in A
                reverse_B_align += self._seqB[j-1]
                j -= 1  # move left
            
            else: # we are in _gapB_matrix
                if self._back_B[i,j] == 0: # move to _align_matrix since that is where we came from when we opened the gap in B
                    ε = 'M'  # move back to _align_matrix
                reverse_A_align += self._seqA[i-1]
                reverse_B_align += '_' # gap in B
                i -= 1 # move up

        # reverse the alignments since we built them in reverse order
        self.seqA_align = reverse_A_align[::-1]
        self.seqB_align = reverse_B_align[::-1]

        return (self.alignment_score, self.seqA_align, self.seqB_align)
    
    def align_matrix(self):
        """
        This function prints the alignment matrix. Useful for debugging.
        """
        print("Alignment Matrix:")
        print(self._align_matrix) 
        return self._align_matrix

    def gapA_matrix(self):
        """
        This function prints the gapA matrix. Useful for debugging.
        """
        print("Gap A Matrix:")
        print(self._gapA_matrix)  
        return self._gapA_matrix

    def gapB_matrix(self):
        """
        This function prints the gapB matrix. Useful for debugging.
        """
        print("Gap B Matrix:")
        print(self._gapB_matrix)
        return self._gapB_matrix


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header


    