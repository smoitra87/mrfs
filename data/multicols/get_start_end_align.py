import os, sys
from Bio import SeqIO

def return_start_end(seq):
    """ Returns the start and end indices of amino acids  """    
    start = next(idx for (idx,s) in enumerate(seq) if s != '-')
    try:
        end = start -1 + next(idx for (idx,s) in enumerate(seq[start:]) if s == '-')
    except StopIteration:
        end = len(seq) - 1

    return start,end

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--alignedf', type=str, help="Aligned Fasta file")
    args = parser.parse_args()
    
    seq_gen = SeqIO.parse(args.alignedf,"fasta")
    baker_seq = str(next(seq_gen).seq)
    pdb_seq = str(next(seq_gen).seq)

    baker_start, baker_end = return_start_end(baker_seq) 
    pdb_start, pdb_end = return_start_end(pdb_seq) 
	
    print baker_start, baker_end
    print pdb_start, pdb_end


