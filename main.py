# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    # VScode suggested most of the below code
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    gg_score, gg_seqA, gg_seqB = nw.align(hs_seq, gg_seq)
    mm_score, mm_seqA, mm_seqB = nw.align(hs_seq, mm_seq)
    br_score, br_seqA, br_seqB = nw.align(hs_seq, br_seq)
    tt_score, tt_seqA, tt_seqB = nw.align(hs_seq, tt_seq)

    # Create a list of tuples (score, species_header) and sort in descending order
    scores = [(gg_score, gg_header), (mm_score, mm_header), (br_score, br_header), (tt_score, tt_header)]
    scores.sort(key=lambda x: x[0], reverse=True) # some fancy sort stuff from VScode

    # Print species in order of most similar to human BRD
    for score, header in scores:
        print(header, '\n', f"Score: {score}\n")

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    # I guess I'll just print the scores again here
    print(f"Human vs Tursiops_truncatus: {tt_score}")
    print(f"Human vs Mus_musculus: {mm_score}")
    print(f"Human vs Gallus_gallus: {gg_score}")
    print(f"Human vs Balaeniceps_rex: {br_score}")
    
    

if __name__ == "__main__":
    main()
